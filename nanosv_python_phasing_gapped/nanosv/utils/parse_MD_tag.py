#!/usr/bin/python
import pysam
import re
import sys
import time
import os

bam = pysam.AlignmentFile(sys.argv[1], 'rb')

for line in bam:
    if line.has_tag('MD'):
        md_tag = re.split("\^\D+",line.get_tag('MD'))
        read_cursor = 0
        ref_cursor = 0
        md_cursor = 0
        md_string = ''
        for t in line.cigartuples:
            if t[0] == 4:
                read_cursor += t[1]
            if t[0] == 1:
                read_cursor += t[1]                
            if t[0] == 2:
                ref_cursor += t[1]
                md_cursor += 1
                md_string = ''               
            if t[0] == 0:
                if re.search("\D+",md_tag[md_cursor]):
                    if md_string == '':
                        for m in re.split("\D",md_tag[md_cursor]):
                            if md_string != '':
                                md_string += 'X'
                            md_string += "="*int(m)                        
                    for seq_match in re.finditer('X',md_string[:t[1]]):
                        print( 'OUT', line.reference_name, (line.reference_start+1+ref_cursor+seq_match.start() ), line.seq[read_cursor+seq_match.start()] )
                    md_string = md_string[t[1]:]
                ref_cursor += t[1]
                read_cursor += t[1]

