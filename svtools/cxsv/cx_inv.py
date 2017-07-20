#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Classification of complex inversion events.
"""

from collections import defaultdict
import numpy as np
import svtools.utils as svu


def breakpoint_ordering(FF, RR, mh_buffer=50):
    """
    Match paired breakpoints to known coordinate orderings.

    e.g. in a simple inversion, FF_start < RR_start < FF_end < RR_end.

    Parameters
    ----------
    FF : pysam.VariantRecord
    RR : pysam.VariantRecord
    mh_buffer : int, optional
        Microhomology buffer. Add window to coordinates to permit fuzzy match.

    Returns
    -------
    ordering : str
        SIMPLE/DEL, DUP5/INS3, DUP3/INS5, dupINVdup, or UNK.
    """

    # Check if breakpoints match simple/deletion ordering
    # (FF_start < RR_start < FF_end < RR_end)
    del_order = ((RR.pos > FF.pos - mh_buffer) and
                 (FF.info['END'] > RR.pos) and
                 (RR.info['END'] > FF.info['END'] - mh_buffer))

    # Check if breakpoints match 5' dup ordering
    # (RR_start < FF_start < FF_end < RR_end)
    dup5_order = ((RR.pos < FF.pos) and
                  (FF.pos < FF.info['END']) and
                  (FF.info['END'] < RR.info['END'] + mh_buffer))

    # Check if breakpoints match 3' dup ordering
    # (FF_start < RR_start < RR_end < FF_end)
    dup3_order = ((FF.pos < RR.pos + mh_buffer) and
                  (RR.pos < RR.info['END']) and
                  (RR.info['END'] < FF.info['END']))

    # Check if breakpoints match dupINVdup ordering
    # (RR_start < FF_start < RR_end < FF_end)
    dupINVdup_order = (RR.pos < FF.pos < RR.info['END'] < FF.info['END'])

    if del_order:
        return 'SIMPLE/DEL'
    elif dup5_order:
        return 'DUP5/INS3'
    elif dup3_order:
        return 'DUP3/INS5'
    elif dupINVdup_order:
        return 'dupINVdup'
    else:
        return 'UNK'


def breakpoints_match(FF, RR, svtype, mh_buffer=50):
    """
    Test if the coordinate ordering matches the class predicted by CNV overlap.

    Parameters
    ----------
    FF : pysam.VariantRecord
    RR : pysam.VariantRecord
    svtype : str
        Class of complex SV predicted by CNV overlap.
        (delINV, INVdel, delINVdel, dupINV, dupINVdel, INVdup, delINVdup,
         dupINVdup)
    mh_buffer : int, optional
        Microhomology buffer. Add window to coordinates to permit fuzzy match.

    Returns
    -------
    ordering : str
        SIMPLE/DEL, DUP5/INS3, DUP3/INS5, dupINVdup, or UNK
    """
    order = breakpoint_ordering(FF, RR, mh_buffer)

    if svtype in 'delINV INVdel delINVdel'.split():
        return order == 'SIMPLE/DEL'
    elif svtype in 'dupINV dupINVdel'.split():
        return order == 'DUP5/INS3'
    elif svtype in 'INVdup delINVdup'.split():
        return order == 'DUP3/INS5'
    else:
        return order == 'dupINVdup'


def classify_2_cnv(FF, RR, cnvs, min_frac=0.5):
    """
    Classify the cxSV class of a pair of inv bkpts and two associated CNVs.

    Matches each CNV to a 5' or 3' location, as constrained by the breakpoint
    coordinates.

    Parameters
    ----------
    FF : pysam.VariantRecord
    RR : pysam.VariantRecord
    cnvs : [pysam.VariantRecord, pysam.VariantRecord]
    min_frac : float, optional
        Minimum reciprocal overlap of each cnv with a candidate CNV interval
        defined by the breakpoint coordinates.

    Returns
    -------
    svtype : str
    """

    # Assign CNVs to 5' or 3' based on ordering
    cnv5, cnv3 = sorted(cnvs, key=lambda r: r.pos)

    # Check if 5' CNV matches breakpoints
    if cnv5.info['SVTYPE'] == 'DEL':
        interval5 = (FF.pos, RR.pos)
    else:
        interval5 = (RR.pos, FF.pos)
    frac5 = svu.reciprocal_overlap(cnv5.pos, cnv5.info['END'], *interval5)

    # Check if 3' CNV matches breakpoints
    if cnv3.info['SVTYPE'] == 'DEL':
        interval3 = (FF.info['END'], RR.info['END'])
    else:
        interval3 = (RR.info['END'], FF.info['END'])
    frac3 = svu.reciprocal_overlap(cnv3.pos, cnv3.info['END'], *interval3)

    # Report cxSV class based on whether CNVs matched intervals
    if frac5 >= min_frac and frac3 >= min_frac:
        svtype = (cnv5.info['SVTYPE'].lower() +
                  'INV' +
                  cnv3.info['SVTYPE'].lower())
    elif frac5 >= min_frac and frac3 < min_frac:
        svtype = classify_1_cnv(FF, RR, cnv5)
    elif frac5 < min_frac and frac3 >= min_frac:
        svtype = classify_1_cnv(FF, RR, cnv3)
    else:
        svtype = 'CNV_2_FAIL'

    return svtype


def classify_1_cnv(FF, RR, cnv, min_frac=0.5,
                   min_bkpt_cnv_size=500, max_bkpt_cnv_size=4000):
    """
    Classify the cxSV class of a pair of inv bkpts and one associated CNV.

    Matches each CNV to a 5' or 3' location, as constrained by the breakpoint
    coordinates. After matching CNV, check if distance between breakpoints
    at other end is sufficient to call a second flanking CNV.

    Parameters
    ----------
    FF : pysam.VariantRecord
    RR : pysam.VariantRecord
    cnvs : [pysam.VariantRecord, pysam.VariantRecord]
    min_frac : float, optional
        Minimum reciprocal overlap of each cnv with a candidate CNV interval
        defined by the breakpoint coordinates.
    min_bkpt_cnv_size : int, optional
        Minimum distance between breakpoints to call flanking CNV.
    max_bkpt_cnv_size : int, optional
        Maximum distance between breakpoints to call flanking CNV.

    Returns
    -------
    svtype : str
    """

    # Make CNV class lowercase (for later concatenation with INV)
    cnv_type = cnv.info['SVTYPE'].lower()

    # Determine eligible 5'/3' CNV intervals defined by the breakpoints
    if cnv_type == 'del':
        interval5 = (FF.pos, RR.pos)
        interval3 = (FF.info['END'], RR.info['END'])
    else:
        interval5 = (RR.pos, FF.pos)
        interval3 = (RR.info['END'], FF.info['END'])

    # Check overlap of CNV against full inversion length
    start = min(FF.pos, RR.pos)
    end = max(FF.info['END'], RR.info['END'])
    total_frac = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], start, end)
    frac5 = svu.overlap_frac(*interval5, cnv.pos, cnv.info['END'])
    frac3 = svu.overlap_frac(*interval3, cnv.pos, cnv.info['END'])

    # If one CNV spans the entire event, it likely represents two CNV merged
    # during preprocessing or clustering
    if total_frac > 0.9 and frac5 > 0.95 and frac3 > 0.95:
        svtype = cnv_type + 'INV' + cnv_type  # + '_merged'
        return svtype

    # Otherwise, check whether it's 5' or 3'
    frac5 = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], *interval5)
    frac3 = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], *interval3)

    # 5' CNV; check 3' breakpoints for small flanking CNV
    if frac5 >= min_frac and frac3 < min_frac:
        svtype = cnv_type + 'INV'

        dist3 = RR.info['END'] - FF.info['END']
        if min_bkpt_cnv_size <= dist3 < max_bkpt_cnv_size:
            svtype = svtype + 'del'
        elif min_bkpt_cnv_size <= -dist3 < max_bkpt_cnv_size:
            svtype = svtype + 'dup'

    # 3' CNV; check 5' breakpoints for small flanking CNV
    elif frac5 < min_frac and frac3 >= min_frac:
        svtype = 'INV' + cnv_type

        dist5 = RR.pos - FF.pos
        if min_bkpt_cnv_size <= dist5 < max_bkpt_cnv_size:
            svtype = 'del' + svtype
        elif min_bkpt_cnv_size <= -dist5 < max_bkpt_cnv_size:
            svtype = 'dup' + svtype

    # Couldn't match the CNV
    else:
        return 'CNV_1_unclassified'

    return svtype


def filter_multiple_cnvs(FF, RR, cnvs, min_frac=0.5):
    """
    For cases with 3 or more overlapping CNV, try to remove spurious hits
    by forcing 50% reciprocal with any of the possible CNV intervals. If
    multiple CNVs are present for a candidate interval (e.g. 5' deletion),
    their coordinates are merged by taking the median.

    Parameters
    ----------
    FF : pysam.VariantRecord
        FF inversion breakpoint.
    RR : pysam.VariantRecord.
        RR inversion breakpoint
    cnvs : list of pysam.VariantRecord
        List of CNVs overlapping breakpoints.

    Returns
    -------
    cnvs : list of pysam.VariantRecord
        Filtered and merged CNVs.
    """

    # Identify eligible intervals for flanking CNV, defined by inv breakpoints
    del5 = (FF.pos, RR.pos)
    del3 = (FF.info['END'], RR.info['END'])
    dup5 = (RR.pos, FF.pos)
    dup3 = (RR.info['END'], FF.info['END'])

    # Determine if CNV supports 5' CNV, 3' CNV, spans event, or fails overlap
    def _test_overlap(cnv):
        svtype = cnv.info['SVTYPE']
        if svtype == 'DEL':
            frac5 = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], *del5)
            frac3 = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], *del3)
        else:
            frac5 = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], *dup5)
            frac3 = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], *dup3)

        if frac5 >= min_frac and frac3 >= min_frac:
            return svtype + '_53'
        elif frac5 >= min_frac:
            return svtype + '_5'
        elif frac3 >= min_frac:
            return svtype + '_3'
        else:
            return 'no_hit'

    # Collect CNV of same overlap type (e.g., 5' deletion) for merging
    cnvlists = defaultdict(list)
    for cnv in cnvs:
        cnvtype = _test_overlap(cnv)
        if cnvtype == 'no_hit':
            continue
        cnvlists[cnvtype].append(cnv)

    # Keep original CNV if only one present,
    # else merge by taking median start/end
    cnvs = []
    for overlap in cnvlists.keys():
        if len(cnvlists[overlap]) == 1:
            cnvs.append(cnvlists[overlap][0])
        else:
            cnvlist = cnvlists[overlap]
            # Overwrite values in first VariantRecord
            # (can't add list of IDs yet)
            merged_cnv = cnvlist[0]

            # get coordinates
            start = int(np.median([c.pos for c in cnvlist]))
            end = int(np.median([c.info['END'] for c in cnvlist]))
            name = '__'.join([c.name for c in cnvlist])

            merged_cnv.pos = start
            merged_cnv.info['END'] = end
            merged_cnv.id = name

            cnvs.append(merged_cnv)

    return sorted(cnvs, key=lambda record: record.pos)


def classify_0_cnv(FF, RR, min_bkpt_cnv_size=300):
    """
    Classify the cxSV class of a pair of inv bkpts with no associated CNVs.

    Matches breakpoint ordering to a known class, then tests if breakpoint
    distances are sufficient to call a CNV that was missed by integration
    pipeline.

    Parameters
    ----------
    FF : pysam.VariantRecord
    RR : pysam.VariantRecord
    min_bkpt_cnv_size : int, optional
        Minimum distance between breakpoints to call flanking CNV.

    Returns
    -------
    svtype : str
    """

    # Identify breakpoint ordering
    order = breakpoint_ordering(FF, RR, mh_buffer=50)

    # Check for flanking deletions around a "simple" inversion
    if order == 'SIMPLE/DEL':
        start_dist = RR.pos - FF.pos
        end_dist = RR.info['END'] - FF.info['END']

        if start_dist < min_bkpt_cnv_size and end_dist < min_bkpt_cnv_size:
            return 'SIMPLE'
        elif start_dist >= min_bkpt_cnv_size and end_dist < min_bkpt_cnv_size:
            return 'delINV'
        elif start_dist < min_bkpt_cnv_size and end_dist >= min_bkpt_cnv_size:
            return 'INVdel'
        else:
            return 'delINVdel'

    # Check for flanking dups
    elif order == 'dupINVdup':
        start_dist = FF.pos - RR.pos
        end_dist = FF.info['END'] - RR.info['END']

        if start_dist >= min_bkpt_cnv_size and end_dist >= min_bkpt_cnv_size:
            return 'dupINVdup'
        elif start_dist >= min_bkpt_cnv_size:
            return 'DUP5/INS3'
        elif end_dist >= min_bkpt_cnv_size:
            return 'DUP3/INS5'
        else:
            return 'UNK'

    # DUP5/INS3, DUP3/INS5, and UNK don't require add'l check
    else:
        return order


def classify_complex_inversion(FF, RR, cnvs):
    """
    Classify the complex class of an inversion and associated CNVs.

    Parameters
    ----------
    FF : pysam.VariantRecord
        FF inversion breakpoint.
    RR : pysam.VariantRecord
        RR inversion breakpoint.
    cnvs : list of pysam.VariantRecord
        List of overlapping CNVs.

    Returns
    -------
    svtype : str
        Complex SV class.
    """

    if len(cnvs) > 2:
        cnvs = filter_multiple_cnvs(FF, RR, cnvs)

    if len(cnvs) == 0:
        #  return 'SIMPLE'
        return classify_0_cnv(FF, RR)
    elif len(cnvs) == 1:
        svtype = classify_1_cnv(FF, RR, cnvs[0])
    elif len(cnvs) == 2:
        svtype = classify_2_cnv(FF, RR, cnvs)
    else:
        return 'MULT_CNVS'

    if breakpoints_match(FF, RR, svtype, mh_buffer=50):
        return svtype
    else:
        return 'COMPLEX_INS'
