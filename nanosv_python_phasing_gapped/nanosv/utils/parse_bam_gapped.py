#!/usr/bin/python
import pysam
import re
import sys
import time

coverages = []
reads = {}
segments = {}
segmentID = 1

bamfile = sys.argv[1]
min_indel_size = '20'

def create_pattern(min_indel_size):    
    pattern_list = [''] * len(min_indel_size)
    for i in range(0,len(min_indel_size)):
        j = len(min_indel_size)-i-1        
        if j > 0:
            pattern_list[i] += '['+str(int(min_indel_size[i])+1)+'-9]\d{'+str(j)+',}'
        for k in range(i+1,len(pattern_list)):
            pattern_list[k] += '['+str(int(min_indel_size[i]))+'-9]'
    pattern_list[-1] += '['+str(int(min_indel_size[-1]))+'-9]'
    pattern = re.compile(r''+"("+"|".join(pattern_list)+")(D)"+'')
    
    return pattern

def search_for_indels(line, clip, clip_2):
    b = False    
    cigarstring = line.cigarstring
    reference_start = line.reference_start+1
    line_query_alignment_length = line.query_alignment_length
    query_alignment_length = 0
    c_start = 0
    pattern = create_pattern(min_indel_size) 
    i = 0
    for m in re.finditer(pattern, line.cigarstring):        
        line.cigarstring = cigarstring[c_start:m.start(0)]
        reference_end = reference_start+line.get_cigar_stats()[0][0]+line.get_cigar_stats()[0][2]+line.get_cigar_stats()[0][7]+line.get_cigar_stats()[0][8]-1
        if line.flag & 16:
            clip_2 = clip_2+query_alignment_length
            query_alignment_length = line.get_cigar_stats()[0][0]+line.get_cigar_stats()[0][1]+line.get_cigar_stats()[0][7]+line.get_cigar_stats()[0][8]            
            if i == 0:
                clip = clip+line_query_alignment_length-query_alignment_length
            else:
                clip = clip-query_alignment_length
            i += 1
        else:
            clip = clip+query_alignment_length
            query_alignment_length = line.get_cigar_stats()[0][0]+line.get_cigar_stats()[0][1]+line.get_cigar_stats()[0][7]+line.get_cigar_stats()[0][8]            
            if i == 0:
                clip_2 = clip_2+line_query_alignment_length-query_alignment_length
            else:
                clip_2 = clip_2-query_alignment_length
            i += 1
            
        print(line.reference_name, reference_start, reference_end, query_alignment_length, clip, clip_2, int(query_alignment_length)+int(clip)+int(clip_2) )
        
        reference_start = int(reference_end)+1+int(m.group(1))
        c_start = m.end(0)
        b = True
    line.cigarstring = cigarstring[c_start:]
    reference_end = reference_start+line.get_cigar_stats()[0][0]+line.get_cigar_stats()[0][2]+line.get_cigar_stats()[0][7]+line.get_cigar_stats()[0][8]-1
    if line.flag & 16:
        clip_2 = clip_2+query_alignment_length
        query_alignment_length = line.get_cigar_stats()[0][0]+line.get_cigar_stats()[0][1]+line.get_cigar_stats()[0][7]+line.get_cigar_stats()[0][8]                    
        clip = clip-query_alignment_length
    else:
        clip = clip+query_alignment_length
        query_alignment_length = line.get_cigar_stats()[0][0]+line.get_cigar_stats()[0][1]+line.get_cigar_stats()[0][7]+line.get_cigar_stats()[0][8]                    
        clip_2 = clip_2-query_alignment_length
    print(line.reference_name, reference_start, reference_end, query_alignment_length, clip, clip_2, int(query_alignment_length)+int(clip)+int(clip_2) )
    
    return b

sys.stderr.write(time.strftime("%c") + " Busy with parsing bam file...\n")
bam = pysam.AlignmentFile(bamfile, 'rb')
if not bam.has_index():
    sys.exit('The bam has no index file')
header = bam.header
if 'HD' in header:
    if not header['HD']['SO'] == 'coordinate':
        sys.exit('The bam file is not coordinate sorted')
if 'RG' in header:
    if type(header['RG']) is list:
        sample_name = header['RG'][0]['SM']
    else:
        sample_name = header['RG']['SM']
else:
    sample_name = re.sub('(\.sorted)?\.bam$', '', str(bamfile))

for line in bam:
    if line.flag & 4 or line.mapping_quality < 20:
        continue
    if line.flag & 16:
        continue
        if line.cigartuples[-1][0] == 5 or line.cigartuples[-1][0] == 4:
            clip = line.cigartuples[-1][1]
        else:
            clip = 0
        if line.cigartuples[0][0] == 5 or line.cigartuples[0][0] == 4:
            clip_2 = line.cigartuples[0][1]
        else:
            clip_2 = 0
    else:
        if line.cigartuples[0][0] == 5 or line.cigartuples[0][0] == 4:
            clip = line.cigartuples[0][1]
        else:
            clip = 0
        if line.cigartuples[-1][0] == 5 or line.cigartuples[-1][0] == 4:
            clip_2 = line.cigartuples[-1][1]
        else:
            clip_2 = 0
    print( clip, clip_2 )
    print( line.query_alignment_length, 3951, 24 )
    b = search_for_indels(line, clip, clip_2)
    
    if b:

        break