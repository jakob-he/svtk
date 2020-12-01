#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Collect split read and discordant pair data from a bam alignment.

Split reads: The tool counts the number of reads soft-clipped in each direction
(30S121M = left-clipped, 121M30S = right-clipped) at each position in the
genome.  The position of a right-clipped read is shifted by the length of its
alignment.

Discordant pairs: The tool reduces discordant pairs to (chrA, posA, strandA,
chrB, posB, strandB).

Unmapped reads, secondary and supplementary
alignments, and duplicates are excluded (SAM flag 256).

Collection can be performed on an S3-hosted bam. The tool will attempt to find
a local copy of the bam index in the working directory, or the directory
specified with `--index-dir`, otherwise the index will be downloaded.
"""

import argparse
import sys
from collections import defaultdict, deque
import numpy as np
import pysam
from natsort import natsorted
import svtk.utils as svu


class SRCollection:
    def __init__(self, bam, splitfile, sample='.',
                 max_split_dist=300):
        self.bam = bam
        self.splitfile = splitfile
        self.discfile = discfile
        self.sample = sample

        # SR evidence
        self.right_split_counts = defaultdict(int)
        self.left_split_counts = defaultdict(int)
        self.prev_split_pos = None
        self.curr_chrom = None
        self.max_split_dist = max_split_dist

    def collect_sr(self):
        """
        Collect SR evidence from a BAM file.

        Reads are considered splitif their CIGAR string contains a soft clip operation.
        """

        for read in self.bam:
            # Restrict to unique primary alignments
            # Equivalent to `samtools view -F 256`
            if svu.is_excluded(read):
                continue

            # Soft clip indicate a candidate split read
            if svu.is_soft_clipped(read):
                if self.splitfile is not None:
                    self.count_split(read)

        self.flush_split_counts()

    def count_split(self, read):
        """
        Count splits at each position.

        Parameters
        ----------
        read : pysam.AlignedSegment
        """

        split_positions = get_split_positions(read)
        #pos, side = get_split_positions(read)

        for (pos, side) in split_positions:
            # Calculate distance to previous split and update position tracker
            # Use abs to catch contig switches
            if self.prev_split_pos is None:
                dist = 0
            else:
                dist = np.abs(pos - self.prev_split_pos)
            self.prev_split_pos = pos

            if self.curr_chrom is None:
                self.curr_chrom = read.reference_name

            # Flush aggregated split reads if we've moved beyond the max dist
            if dist > self.max_split_dist:
                self.flush_split_counts()
                self.curr_chrom = read.reference_name

            # Tally the split at its corresponding position
            if side == 'RIGHT':
                self.right_split_counts[pos] += 1
            elif side == 'LEFT':
                self.left_split_counts[pos] += 1

    def flush_split_counts(self):
        """
        Write current split counts to disk and reset dictionaries
        """

        # Compile counts collected so far
        entries = deque()
        for clip in 'left right'.split():
            df = getattr(self, '%s_split_counts' % clip)

            for pos, count in df.items():
                entries.append((self.curr_chrom, pos, clip, count,
                                self.sample))

        # Sort in chunks as we go
        entries = sorted(entries, key=lambda s: s[1])

        # Flush to disk
        fmt = '%s\t%d\t%s\t%d\t%s\n'
        for entry in entries:
            self.splitfile.write((fmt % entry).encode('utf-8'))

        # Reset split counts
        self.right_split_counts = defaultdict(int)
        self.left_split_counts = defaultdict(int)


def get_split_positions(read):
    """
    Calculate split coordinate based on read alignment and CIGAR operations.

    Support is only present for reads soft-clipped on one side, e.g. 100M51S,
    as the coordinate is calculated by shifting the alignment position by the
    length of the flanking match operation.

    Parameters
    ----------
    read : pysam.AlignedSegment

    Returns
    -------
    pos : int
        Adjusted split read coordinate
    side : str [RIGHT,LEFT,MIDDLE]
        Direction of soft clip
    """


    pos = read.pos

    split_positions = []

    # Left soft clip - sequence is already aligned to split position
    if is_left_clipped(read):
        split_positions.append([pos, 'LEFT'])

    # Right soft clip - add length of aligned sequence
    if is_right_clipped(read):
        clip_pos = pos
        for operation, length in read.cigartuples:
            # Only shift based on matches, ignore DEL/INS/clips
            if not is_clipping_operation(operation) and operation_consumes_ref_bases(operation):
                clip_pos += length
        split_positions.append([clip_pos, 'RIGHT'])

    return split_positions

def is_left_clipped(read):
    return len(read.cigartuples) >= 1 and is_clipping_operation(read.cigartuples[0][0])

def is_right_clipped(read):
    return len(read.cigartuples) >= 1 and is_clipping_operation(read.cigartuples[-1][0])

def is_clipping_operation(operation):
    return operation == 4 or operation == 5

def operation_consumes_ref_bases(operation):
    """
    Returns true if this is a cigar operation that consumes reference bases
    """
    return operation == 0 or operation == 2 or operation == 3 or operation == 7

def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk collect-sr',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('bam', help='Local or S3 path to bam')
    parser.add_argument('sample', help='ID to append to each line of output '
                        'files.')
    parser.add_argument('splitfile',
                        help='Output split counts.')

    parser.add_argument('--index-dir', default=None,
                        help='Directory of local BAM indexes if accessing '
                        'a remote S3 bam.')
    parser.add_argument('-r', '--region',
                        help='Tabix-formatted region to parse')
    parser.add_argument('-z', '--bgzip', default=False, action='store_true',
                        help='bgzip and tabix index output')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Load bam from S3 if necessary
    if args.bam.startswith('s3://'):
        bam = svu.load_s3bam(args.bam, args.index_dir)
    else:
        bam = pysam.AlignmentFile(args.bam)

    # Restrict to region of interest
    if args.region:
        bam = bam.fetch(args.region.encode('utf-8'))

    # Collect data and save
    with svu.BgzipFile(args.splitfile, args.bgzip) as splitfile:
        with svu.BgzipFile(args.discfile, args.bgzip) as discfile:
            SRCollection(bam, splitfile, discfile, args.sample).collect_sr()


if __name__ == '__main__':
    main(sys.argv[1:])
