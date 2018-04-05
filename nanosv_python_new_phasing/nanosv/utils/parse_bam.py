#!/usr/bin/python
import pysam
import re
import sys
import time
import os

from classes import read as r
from classes import segment as s

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NanoSV

coverages = []
reads = {}
segments = {}
segmentID = 1

def parse_bam():
    """
    Reads bam file and saves reads and their segments in objects of the Read en Segment classes.
    :param bamfile used to open bam file:
    """
    global sample_name, header, segmentID, bam, gevonden_varianten, aantal
    aantal = 0
    print("Parse_bam.py gestart")
    sys.stderr.write(time.strftime("%c") + " Busy with parsing bam file...\n")
    bam = pysam.AlignmentFile(NanoSV.opts_bam, 'rb')
    header = bam.header
    if not bam.has_index():
        sys.exit('The bam has no index file')
    if 'HD' in header:
        if not header['HD']['SO'] == 'coordinate':
            sys.exit('The bam file is not coordinate sorted')
    if 'RG' in header:
        sample_name = header['RG']['SM']
    else:
        sample_name = re.sub('(\.sorted)?\.bam$', '', str(NanoSV.opts_bam))
    previous_cursor = -1
    segments_to_check = {}
    gevonden_varianten = []
    for line in bam:
        remove = [end for end in segments_to_check if line.reference_start > end]
        for segment_ID in remove: del segments_to_check[segment_ID]
        if line.query_name in reads:
            read = reads[line.query_name]
        else:
            read = r.Read(line.query_name, line.infer_read_length())
            reads[line.query_name] = read
        
        if line.flag & 4 or line.mapping_quality < NanoSV.opts_min_mapq:
            continue
        segment = s.Segment(segmentID, line.query_name, line.flag, line.reference_name, line.reference_start+1, line.mapping_quality,
                            line.query_alignment_length, find_variations(line.cigartuples, line.pos, line.query_qualities, line.seq))
        segment.end = line.reference_start + line.reference_length
        segment.pid = format(line.get_cigar_stats()[0][7] / segment.length, '.3f')
        if previous_cursor == -1:
            previous_cursor = line.reference_start
        if segment.pid == "0.000":
            segment.pid = format(line.get_cigar_stats()[0][0] / segment.length, '.3f')
        if line.flag & 16:
            if line.cigartuples[-1][0] == 5 or line.cigartuples[-1][0] == 4:
                segment.clip = line.cigartuples[-1][1]
            else:
                segment.clip = 0
            if line.cigartuples[0][0] == 5 or line.cigartuples[0][0] == 4:
                segment.clip_2 = line.cigartuples[0][1]
            else:
                segment.clip_2 = 0
        else:
            if line.cigartuples[0][0] == 5 or line.cigartuples[0][0] == 4:
                segment.clip = line.cigartuples[0][1]
            else:
                segment.clip = 0
            if line.cigartuples[-1][0] == 5 or line.cigartuples[-1][0] == 4:
                segment.clip_2 = line.cigartuples[-1][1]
            else:
                segment.clip_2 = 0
        if float(segment.pid) < NanoSV.opts_min_pid:
            continue
        read.addSegment(segment)
        segments[segmentID] = segment
        # if line.reference_start > previous_cursor:
        #     remove_variations(previous_cursor, line.reference_start, segments_to_check)
        segments_to_check[segmentID] = segment.end
        segmentID += 1
    # print(len(gevonden_varianten), aantal)


def find_variations(cigartuples, pos, qual, seq):
    global gevonden_varianten, aantal
    cursor = (int(pos) - 1)
    seq_base_index = 0
    variations = {}
    base_to_number = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    for tuple in cigartuples:
        if tuple[0] == 8:
            for mismatch in range(tuple[1]):
                cursor += 1
                variations[cursor] = [base_to_number[seq[seq_base_index]], qual[seq_base_index]]
                seq_base_index += 1
                gevonden_varianten.append(cursor)
                aantal += 1
        elif tuple[0] == 2:
            for i in range(tuple[1]):
                cursor += 1
                variations[cursor] = [-1]
        elif tuple[0] == 7:
            for equal in range(tuple[1]):
                cursor += 1
                seq_base_index += 1
        elif tuple[0] == 1:
            for ins in range(tuple[1]):
                seq_base_index += 1
    return variations


def remove_variations(start, end, segments_to_check):
    global gevonden_varianten
    for pos in range(start, end):
        values_of_position_highq = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
        values_of_position_lowq = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
        for segment in segments_to_check:
            # print(segments[segment].variations)
            try:
                if segments[segment].variations[pos] == [0]:
                    values_of_position_highq[0] += 1
                elif segments[segment].variations[pos][1] >= NanoSV.opts_min_base_qual_ph:
                    values_of_position_highq[segments[segment].variations[pos][0]] += 1
                elif segments[segment].variations[pos][1] < NanoSV.opts_min_base_qual_ph:
                    values_of_position_lowq[segments[segment].variations[pos][0]] += 1
            except KeyError:
                values_of_position_highq[5] += 1
        valuable_call = 0
        if values_of_position_highq[0] > (NanoSV.opts_max_deletions * len(segments_to_check)):
            valuable_call -= 5
        for key, value in values_of_position_highq.items():
            if key == 0:
                continue
            if key != 5 and value > NanoSV.opts_min_occurences_of_var * len(segments_to_check) and value + \
                    values_of_position_lowq[key] < int(NanoSV.opts_max_occurences_of_var * len(segments_to_check)):
                valuable_call += 1
        if valuable_call < 1:
            if pos in gevonden_varianten:

                del gevonden_varianten[gevonden_varianten.index(pos)]
            for segment in segments_to_check:
                try:
                    del segments[segment].variations[pos]
                except KeyError:
                    print("keyerror bij delete")
                    continue
