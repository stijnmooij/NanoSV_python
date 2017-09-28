import optparse
import os
import sys
import textwrap
import time
import json

from collections import defaultdict

import read as r
import segment as s
import breakpoint as b
import sv as svclass

import re

parser = optparse.OptionParser()
group1 = optparse.OptionGroup(parser, 'General options')
group2 = optparse.OptionGroup(parser, 'Filter options')
group3 = optparse.OptionGroup(parser, 'Detection options')
group4 = optparse.OptionGroup(parser, 'Output filter options')
group1.add_option('-t', '--threads', dest='threads', help='[i] Number of threads. Default: 8')
group1.add_option('--sambamba', dest='sambamba', help='[s] Path to sambamba: Default: sambamba_v0.6.3')
group2.add_option('-s', '--split', dest='split', help='[i]Maximum number of segments per read. Default: 10')
group2.add_option('-p', '--pid', dest='pid', help='[s] Minimum percentage identity to reference. Default: 0.70')
group2.add_option('-m', '--mapq', dest='mapq', help='[i] Minimum mapping quality. Default: 20')
group3.add_option('-d', '--distance', dest='distance', help='[i] Maximum distance to cluster SVs together. Default: 10')
group3.add_option('-c', '--count', dest='count', help='[i] Minimum number of supporting reads. Default: 2')
group3.add_option('-f', '--refdistance', dest='refdistance', help='Minimum distance for reference reads: Default: 100')
group3.add_option('-u', '--unmapped', dest='unmapped', help='[i] Minimum unmapped length of hanging segments: 20')
group3.add_option('-r', '--matedistance', dest='matedistance',
                  help='[i] Maximum distance to look for mateid. Default: 300')
group4.add_option('-w', '--window', dest='window', help='[i] Maximum window size. Default: 1000')
group4.add_option('-n', '--cluster', dest='cluster', help='[i] Maximum number of SVs in a window. Default: 2')
group4.add_option('-q', '--mapqf', dest='mapqf', help='Minimum median mapping quality of a SV. Default: 80')
group4.add_option('-i', '--pidf', dest='pidf',
                  help='[s] Minimum median percentage identity to reference. Default: 0.80')
group4.add_option('-g', '--gap', dest='gap', help='[i] Maximum median gap size. Default: 100')
group4.add_option('-y', '--ci', dest='ci', help='[i] Maximum Confidence interval distance. Default: 20')
parser.add_option_group(group1)
parser.add_option_group(group2)
parser.add_option_group(group3)
parser.add_option_group(group4)
(opts, args) = parser.parse_args()

bam = sys.argv[1]

opts.threads = 8
opts.sambamba = '/data/sambamba_v0.6.3'
opts.split = 10
opts.pid = 0.70
opts.mapq = 20
opts.distance = 10
opts.count = 2
opts.refdistance = 100
opts.unmapped = 20
opts.matedistance = 300
opts.window = 1000
opts.cluster = 2
opts.mapqf = 80
opts.pidf = 0.80
opts.gap = 100
opts.ci = 20

reads = {}
segments = {}
breakpoints = {}
breakpoints_region = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
hanging_breakpoints_region = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
structural_variants = {}
structural_variants_region = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
structural_variants_region_2 = defaultdict(lambda: defaultdict(int))

aantal_svs = 0
svID = 1
segmentID = 1


def parse_svs():
    global aantal_svs
    print(time.strftime("%c") + " Busy with reference reads...")
    prev_rname = 0
    for segment_id in sorted(segments):
        if segments[segment_id].rname != prev_rname and prev_rname != 0:
            del structural_variants_region[prev_rname]
        for sv_min in sorted(structural_variants_region[segments[segment_id].rname]):
            if sv_min <= segments[segment_id].pos:
                del structural_variants_region[segments[segment_id].rname][sv_min]
                continue
            if sv_min > segments[segment_id].end:
                break
            for sv_max in sorted(structural_variants_region[segments[segment_id].rname][sv_min]):
                if sv_max >= segments[segment_id].end:
                    break
                for sv_id in structural_variants_region[segments[segment_id].rname][sv_min][sv_max]:
                    x = structural_variants_region[segments[segment_id].rname][sv_min][sv_max][sv_id]
                    if sv_min > (segments[segment_id].pos + opts.refdistance) and sv_max < (
                                segments[segment_id].end - opts.refdistance):
                        structural_variants[sv_id].format['DR'][x] += 1
                        structural_variants[sv_id].format['RO'][x] += (1 - 10 ** (-segments[segment_id].mapq / 10.0))
        prev_rname = segments[segment_id].rname

    print(time.strftime("%c") + " Busy with printing to vcf...")
    for sv_id in sorted(structural_variants):
        aantal_svs += 1
        if structural_variants[sv_id].set != 1: structural_variants[sv_id].setArguments()
        for sv_id_2 in sorted(structural_variants_region_2[structural_variants[sv_id].chr]):
            structural_variants[sv_id].setCluster(1)
            structural_variants[sv_id_2].setCluster(1)
            if sv_id > sv_id_2: continue
            if structural_variants[sv_id_2].set != 1: structural_variants[sv_id_2].setArguments()
            if structural_variants[sv_id].chr != structural_variants[sv_id_2].chr: continue
            if structural_variants[sv_id].chr2 != structural_variants[sv_id_2].chr2: continue
            if abs(structural_variants[sv_id].pos - structural_variants[sv_id_2].pos) <= opts.window and abs(
                            structural_variants[sv_id].info['END'] - structural_variants[sv_id_2].info[
                        'END']) <= opts.window:
                structural_variants[sv_id].SVcluster += 1
                structural_variants[sv_id_2].SVcluster += 1
            if structural_variants[sv_id_2].info['SVTYPE'] != "BND": continue
            if structural_variants[sv_id].info['SVTYPE'] != "BND": continue
            if re.match("/\](\w+):(\d+)\]\w+/", structural_variants[sv_id].alt) and not re.match("/\w+\[(\w+):(\d+)\[/",
                                                                                                 structural_variants[
                                                                                                     sv_id_2].alt): continue
            if re.match("/\[(\w+):(\d+)\[\w+/", structural_variants[sv_id].alt) and not re.match("/\w+\](\w+):(\d+)\]/",
                                                                                                 structural_variants[
                                                                                                     sv_id_2].alt): continue
            if re.match("/\w+\](\w+):(\d+)\]/", structural_variants[sv_id].alt) and not re.match("/\[(\w+):(\d+)\[\w+/",
                                                                                                 structural_variants[
                                                                                                     sv_id_2].alt): continue
            if re.match("/\w+\[(\w+):(\d+)\[/", structural_variants[sv_id].alt) and not re.match("/\](\w+):(\d+)\]\w+/",
                                                                                                 structural_variants[
                                                                                                     sv_id_2].alt): continue
            if abs(structural_variants[sv_id].pos - structural_variants[sv_id_2].pos) > opts.matedistance: continue
            if abs(structural_variants[sv_id].info['END'] - structural_variants[sv_id_2].info[
                'END']) > opts.matedistance: continue
            structural_variants[sv_id].info['MATEID'] = sv_id_2
            structural_variants[sv_id_2].info['MATEID'] = sv_id
        if structural_variants[sv_id].SVcluster > opts.cluster: structural_variants[sv_id].filter.append("SVcluster")
        if structural_variants[sv_id].info['GAP'] > opts.gap and structural_variants[sv_id].info['SVTYPE'] != "INS":
            structural_variants[sv_id].filter.append("GAP")
        print(type(structural_variants[sv_id].info['MAPQ']))
        if re.match("/(\d+),(\d+)/", structural_variants[sv_id].info['MAPQ']):
            ma = re.search("/(\d+),(\d+)/", )
            print(ma.group(1), opts.mapqf)
            if ma.group(1) < opts.mapqf or ma.group(2) < opts.mapqf: structural_variants[sv_id].filter = "MapQual"
        if re.match("/(\d.\d+),(\d.\d+)/", structural_variants[sv_id].info['PID']):
            ma = re.search("/(\d.\d+),(\d.\d+)/", )
            if ma.group(1) < opts.pidf or ma.group(2) < opts.pidf: structural_variants[sv_id].filter = "PID"
        if re.match("/(\d+),(\d+)/", structural_variants[sv_id].info['CIPOS']):
            ma = re.search("/(\d+),(\d+)/", )
            if ma.group(1) > opts.ci or ma.group(2) > opts.ci: structural_variants[sv_id].filter = "CIPOS"
        if re.match("/(\d+),(\d+)/", structural_variants[sv_id].info['CIEND']):
            ma = re.search("/(\d+),(\d+)/", )
            if ma.group(1) > opts.ci or ma.group(2) > opts.ci: structural_variants[sv_id].filter = "CIEND"

        if structural_variants[sv_id].info['SVTYPE'] == 'INS': structural_variants[sv_id].info['SVLEN'] = \
            structural_variants[sv_id].info['GAP']
        structural_variants[sv_id].printVCF()


def parse_breakpoints_2(breakpoints_region_2):
    prev_pos_1 = -1
    prev_pos_2 = -1
    global svID
    for pos_2 in sorted(breakpoints_region_2):
        for pos_1 in sorted(breakpoints_region_2[pos_2]):
            for breakpoint_id in breakpoints_region_2[pos_2][pos_1]:
                breakpoint = breakpoints[breakpoint_id]
                if (prev_pos_2 == -1):
                    sv = svclass.SV(svID, breakpoint)
                    svID += 1
                elif abs(breakpoint.segment_2["pos"] - prev_pos_2) <= opts.distance:
                    sv.addBreakpoint(breakpoint)
                else:
                    if sum(sv.format['DV']) >= opts.count * 2:
                        if sv.flag1 == 16:
                            for hanging_pos in sorted(hanging_breakpoints_region[sv.chr]['H']):
                                if hanging_pos < (min(sv.pos) - opts.distance):
                                    continue
                                if hanging_pos > (max(sv.pos) + opts.distance):
                                    break
                                for hanging_id in hanging_breakpoints_region[sv.chr]['H'][hanging_pos]:
                                    segment_id = hanging_breakpoints_region[sv.chr]['H'][hanging_pos][hanging_id]
                                    if segment_id == 0: continue
                                    sv.format['HR'][0] += 1
                                    sv.format['VO'][0] += (1 - 10 ** (-segments[segment_id].mapq / 10.0))
                                    sv.addInfoField("PID", [[segments[segment_id].pid], None])
                                    sv.addInfoField("MAPQ", [[segments[segment_id].mapq], None])
                                    sv.addInfoField("PLENGTH", [[segments[segment_id].plength], None])
                                    sv.addInfoField("RLENGTH", [reads[segments[segment_id].qname].length])
                                    if re.match("/2D_2d$/", segments[segment_id].qname):
                                        sv.addInfoField("RT", ["2d"])
                                    elif re.match("/2D_complement$/", segments[segment_id].qname):
                                        sv.addInfoField("RT", ["complement"])
                                    else:
                                        sv.addInfoField("RT", ["template"])
                        else:
                            for hanging_pos in sorted(hanging_breakpoints_region[sv.chr]['T']):
                                if hanging_pos < (min(sv.pos) - opts.distance):
                                    continue
                                if hanging_pos > (max(sv.pos) + opts.distance):
                                    break
                                for hanging_id in hanging_breakpoints_region[sv.chr]['T'][hanging_pos]:
                                    segment_id = hanging_breakpoints_region[sv.chr]['T'][hanging_pos][hanging_id]
                                    if segment_id == 0: continue
                                    sv.format['HR'][0] += 1
                                    sv.format['VO'][0] += (1 - 10 ** (-segments[segment_id].mapq / 10.0))
                                    sv.addInfoField("PID", [[segments[segment_id].pid], None])
                                    sv.addInfoField("MAPQ", [[segments[segment_id].mapq], None])
                                    sv.addInfoField("PLENGTH", [[segments[segment_id].plength], None])
                                    sv.addInfoField("RLENGTH", [reads[segments[segment_id].qname].length])
                                    if re.match("/2D_2d$/", segments[segment_id].qname):
                                        sv.addInfoField("RT", ["2d"])
                                    elif re.match("/2D_complement$/", segments[segment_id].qname):
                                        sv.addInfoField("RT", ["complement"])
                                    else:
                                        sv.addInfoField("RT", ["template"])
                        if sv.flag2 == 16:
                            for hanging_pos in sorted(hanging_breakpoints_region[sv.chr2]['T']):
                                if hanging_pos < (min(sv.info['END']) - opts.distance):
                                    continue
                                if hanging_pos > (max(sv.info['END']) + opts.distance):
                                    break
                                for hanging_id in hanging_breakpoints_region[sv.chr2]['T'][hanging_pos]:
                                    segment_id = hanging_breakpoints_region[sv.chr]['T'][hanging_pos][hanging_id]
                                    if segment_id == 0: continue
                                    sv.format['HR'][1] += 1
                                    sv.format['VO'][0] += (1 - 10 ** (-segments[segment_id].mapq / 10.0))
                                    sv.addInfoField("PID", [None, [segments[segment_id].pid]])
                                    sv.addInfoField("MAPQ", [None, [segments[segment_id].mapq]])
                                    sv.addInfoField("PLENGTH", [None, [segments[segment_id].plength]])
                                    sv.addInfoField("RLENGTH", [reads[segments[segment_id].qname].length])
                                    if re.match("/2D_2d$/", segments[segment_id].qname):
                                        sv.addInfoField("RT", ["2d"])
                                    elif re.match("/2D_complement$/", segments[segment_id].qname):
                                        sv.addInfoField("RT", ["complement"])
                                    else:
                                        sv.addInfoField("RT", ["template"])
                        else:
                            for hanging_pos in sorted(hanging_breakpoints_region[sv.chr2]['H']):
                                if hanging_pos < (min(sv.pos) - opts.distance):
                                    continue
                                if hanging_pos > (max(sv.pos) + opts.distance):
                                    break
                                for hanging_id in hanging_breakpoints_region[sv.chr2]['H'][hanging_pos]:
                                    segment_id = hanging_breakpoints_region[sv.chr]['H'][hanging_pos][hanging_id]
                                    if segment_id == 0: continue
                                    sv.format['HR'][1] += 1
                                    sv.format['VO'][1] += (1 - 10 ** (-segments[segment_id].mapq / 10.0))
                                    sv.addInfoField("PID", [None, [segments[segment_id].pid]])
                                    sv.addInfoField("MAPQ", [None, [segments[segment_id].mapq]])
                                    sv.addInfoField("PLENGTH", [None, [segments[segment_id].plength]])
                                    sv.addInfoField("RLENGTH", [reads[segments[segment_id].qname].length])
                                    if re.match("/2D_2d$/", segments[segment_id].qname):
                                        sv.addInfoField("RT", ["2d"])
                                    elif re.match("/2D_complement$/", segments[segment_id].qname):
                                        sv.addInfoField("RT", ["complement"])
                                    else:
                                        sv.addInfoField("RT", ["template"])
                        structural_variants[sv.id] = sv
                        structural_variants_region[sv.chr][min(sv.pos)][max(sv.pos)][sv.id] = 0
                        structural_variants_region[sv.chr2][min(sv.info['END'])][max(sv.info['END'])][sv.id] = 1
                        structural_variants_region_2[sv.chr][sv.id] = 1

                    sv = svclass.SV(svID, breakpoint)
                    svID += 1

                sv.addInfoField("PID",
                                [[segments[breakpoint.segment_1["id"]].pid], [segments[breakpoint.segment_2["id"]].pid]])
                sv.addInfoField("MAPQ",
                                [[segments[breakpoint.segment_1["id"]].mapq], [segments[breakpoint.segment_2["id"]].mapq]])
                sv.addInfoField("PLENGTH", [[segments[breakpoint.segment_1["id"]].plength],
                                            [segments[breakpoint.segment_2["id"]].plength]])
                sv.addInfoField("RLENGTH", [reads[segments[breakpoint.segment_2["id"]].qname].length])
                sv.addInfoField("GAP", [breakpoint.gap])

                if re.match("/2D_2d$/", segments[breakpoint.segment_2["id"]].qname):
                    sv.addInfoField("RT", ["2d"])
                elif re.match("/2D_complement$/", segments[breakpoint.segment_2["id"]].qname):
                    sv.addInfoField("RT", ["complement"])
                else:
                    sv.addInfoField("RT", ["template"])

                prev_pos_1 = breakpoint.segment_1["pos"]
                prev_pos_2 = int(pos_2)

    if sum(sv.format['DV']) >= opts.count * 2:
        if sv.flag1 == 16:
            for hanging_pos in sorted(hanging_breakpoints_region[sv.chr]['H']):
                if hanging_pos < (min(sv.pos) - opts.distance):
                    continue
                if hanging_pos > (max(sv.pos) + opts.distance):
                    break
                for hanging_id in hanging_breakpoints_region[sv.chr]['H'][hanging_pos]:
                    segment_id = hanging_breakpoints_region[sv.chr]['H'][hanging_pos][hanging_id]
                    if segment_id == 0: continue
                    sv.format['HR'][0] += 1
                    sv.format['VO'][0] += (1 - 10 ** (-segments[segment_id].mapq / 10.0))
                    sv.addInfoField("PID", [[segments[segment_id].pid], None])
                    sv.addInfoField("MAPQ", [[segments[segment_id].mapq], None])
                    sv.addInfoField("PLENGTH", [[segments[segment_id].plength], None])
                    sv.addInfoField("RLENGTH", [reads[segments[segment_id].qname].length])
                    if re.match("/2D_2d$/", segments[segment_id].qname):
                        sv.addInfoField("RT", ["2d"])
                    elif re.match("/2D_complement$/", segments[segment_id].qname):
                        sv.addInfoField("RT", ["complement"])
                    else:
                        sv.addInfoField("RT", ["template"])
        else:
            for hanging_pos in sorted(hanging_breakpoints_region[sv.chr]['T']):
                if hanging_pos < (min(sv.pos) - opts.distance):
                    continue
                if hanging_pos > (max(sv.pos) + opts.distance):
                    break
                for hanging_id in hanging_breakpoints_region[sv.chr]['T'][hanging_pos]:
                    segment_id = hanging_breakpoints_region[sv.chr]['T'][hanging_pos][hanging_id]
                    if segment_id == 0: continue
                    sv.format['HR'][0] += 1
                    sv.format['VO'][0] += (1 - 10 ** (-segments[segment_id].mapq / 10.0))
                    sv.addInfoField("PID", [[segments[segment_id].pid], None])
                    sv.addInfoField("MAPQ", [[segments[segment_id].mapq], None])
                    sv.addInfoField("PLENGTH", [[segments[segment_id].plength], None])
                    sv.addInfoField("RLENGTH", [reads[segments[segment_id].qname].length])
                    if re.match("/2D_2d$/", segments[segment_id].qname):
                        sv.addInfoField("RT", ["2d"])
                    elif re.match("/2D_complement$/", segments[segment_id].qname):
                        sv.addInfoField("RT", ["complement"])
                    else:
                        sv.addInfoField("RT", ["template"])
        if sv.flag2 == 16:
            for hanging_pos in sorted(hanging_breakpoints_region[sv.chr2]['T']):
                if hanging_pos < (min(sv.info['END']) - opts.distance):
                    continue
                if hanging_pos > (max(sv.info['END']) + opts.distance):
                    break
                for hanging_id in hanging_breakpoints_region[sv.chr2]['T'][hanging_pos]:
                    segment_id = hanging_breakpoints_region[sv.chr]['T'][hanging_pos][hanging_id]
                    if segment_id == 0: continue
                    sv.format['HR'][1] += 1
                    sv.format['VO'][1] += (1 - 10 ** (-segments[segment_id].mapq / 10.0))
                    sv.addInfoField("PID", [None, [segments[segment_id].pid]])
                    sv.addInfoField("MAPQ", [None, [segments[segment_id].mapq]])
                    sv.addInfoField("PLENGTH", [None, [segments[segment_id].plength]])
                    sv.addInfoField("RLENGTH", [reads[segments[segment_id].qname].length])
                    if re.match("/2D_2d$/", segments[segment_id].qname):
                        sv.addInfoField("RT", ["2d"])
                    elif re.match("/2D_complement$/", segments[segment_id].qname):
                        sv.addInfoField("RT", ["complement"])
                    else:
                        sv.addInfoField("RT", ["template"])
        else:
            for hanging_pos in sorted(hanging_breakpoints_region[sv.chr2]['H']):
                if hanging_pos < (min(sv.pos) - opts.distance):
                    continue
                if hanging_pos > (max(sv.pos) + opts.distance):
                    break
                for hanging_id in hanging_breakpoints_region[sv.chr2]['H'][hanging_pos]:
                    segment_id = hanging_breakpoints_region[sv.chr]['H'][hanging_pos][hanging_id]
                    if segment_id == 0: continue
                    sv.format['HR'][1] += 1
                    sv.format['VO'][1] += (1 - 10 ** (-segments[segment_id].mapq / 10.0))
                    sv.addInfoField("PID", [None, [segments[segment_id].pid]])
                    sv.addInfoField("MAPQ", [None, [segments[segment_id].mapq]])
                    sv.addInfoField("PLENGTH", [None, [segments[segment_id].plength]])
                    sv.addInfoField("RLENGTH", [reads[segments[segment_id].qname].length])
                    if re.match("/2D_2d$/", segments[segment_id].qname):
                        sv.addInfoField("RT", ["2d"])
                    elif re.match("/2D_complement$/", segments[segment_id].qname):
                        sv.addInfoField("RT", ["complement"])
                    else:
                        sv.addInfoField("RT", ["template"])
        structural_variants[sv.id] = sv
        structural_variants_region[sv.chr][min(sv.pos)][max(sv.pos)][sv.id] = 0
        structural_variants_region[sv.chr2][min(sv.info['END'])][max(sv.info['END'])][sv.id] = 1
        structural_variants_region_2[sv.chr][sv.id] = 1


def parse_breakpoints():
    print(time.strftime("%c") + " Busy with parsing breakpoints...")
    for region in breakpoints_region:
        prev_pos_1 = -1
        breakpoints_region_2 = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        for pos_1 in sorted(breakpoints_region[region]):
            for breakpoint_id in breakpoints_region[region][pos_1]:
                if prev_pos_1 == -1:
                    breakpoints_region_2[breakpoints[breakpoint_id].segment_2["pos"]][pos_1][breakpoint_id] = 1
                elif abs(pos_1 - prev_pos_1) <= opts.distance:
                    breakpoints_region_2[breakpoints[breakpoint_id].segment_2["pos"]][pos_1][breakpoint_id] = 1
                else:
                    parse_breakpoints_2(breakpoints_region_2)
                    breakpoints_region_2.clear()
                    breakpoints_region_2[breakpoints[breakpoint_id].segment_2["pos"]][pos_1][breakpoint_id] = 1
                prev_pos_1 = int(pos_1)
        parse_breakpoints_2(breakpoints_region_2)


def parse_reads():
    print(time.strftime("%c") + " Busy with parsing read segments...")
    breakpointID = 1
    hanging_breakpointID = -1
    for qname in reads:
        clips = sorted(reads[qname].segments)
        if len(clips) == 1 or len(clips) > opts.split:
            continue
        for i in range(0, (len(clips) - 1)):
            i2 = i + 1

            segment_1 = segments[reads[qname].segments[clips[i]]]
            segment_2 = segments[reads[qname].segments[clips[i2]]]
            segment_1.setPlength(reads[qname].length)
            segment_2.setPlength(reads[qname].length)

            breakpoint = b.Breakpoint(breakpointID, segment_1, segment_2)
            breakpointID += 1

            gap = (clips[i2] - (clips[i] + segment_1.length))
            breakpoint.setGap(gap)
            breakpoint.setBreakpoint(segment_1, segment_2)
            if segment_1.rname > segment_2.rname:
                breakpoint.switchSegments()
            elif segment_1.rname == segment_2.rname and breakpoint.segment_1["pos"] > breakpoint.segment_2["pos"]:
                breakpoint.switchSegments()

            breakpoint.setSVtype()
            breakpoints[breakpoint.id] = breakpoint
            values = (breakpoint.svtype, str(breakpoint.segment_1["rname"]), str(breakpoint.segment_2["rname"]),
                      str(breakpoint.segment_1["flag"]), str(breakpoint.segment_2["flag"]))
            breakpoints_region["\t".join(values)][breakpoint.segment_1["pos"]][breakpoint.id] = 1
            if i == 0 and segment_1.clip >= opts.unmapped:
                hanging_breakpoint_pos = segment_1.pos
                if segment_1.flag == 16:
                    hanging_breakpoint_pos = segment_1.end
                    hanging_breakpoints_region[segment_1.rname]['T'][hanging_breakpoint_pos][hanging_breakpointID] = segment_1.id
                    hanging_breakpointID -= 1
                else:
                    hanging_breakpoints_region[segment_1.rname]['H'][hanging_breakpoint_pos][hanging_breakpointID] = segment_1.id
                    hanging_breakpointID = hanging_breakpointID - 1
            if i2 == len(clips) - 1 and segment_2.clip_2 >= opts.unmapped:
                hanging_breakpoint_pos = segment_2.end
                if segment_2.flag == 16:
                    hanging_breakpoint_pos = segment_2.pos
                    hanging_breakpoints_region[segment_2.rname]['H'][hanging_breakpoint_pos][hanging_breakpointID] = segment_2.id
                    hanging_breakpointID -= 1
                else:
                    hanging_breakpoints_region[segment_2.rname]['T'][hanging_breakpoint_pos][hanging_breakpointID] = segment_2.id

                    hanging_breakpointID -= 1


def parse_cigar(cigar):
    cigar_dict = dict()
    for char in list('H=DI'):
        m = re.findall(r"(\d+)" + char, cigar)
        cigar_dict[char] = sum(map(int, m))
    m1 = re.findall(r"^(\d+)H", cigar)
    m2 = re.findall(r"(\d+)H$", cigar)
    cigar_dict['H'] = [sum(map(int, m1)), sum(map(int, m2))]
    return cigar_dict


def parse_bam(bam, segmentID):
    print(time.strftime("%c") + " Busy with parsing bam file...")

    with os.popen(opts.sambamba + ' view -t ' + str(opts.threads) + ' ' + bam) as bam:
        for line in bam:
            line = line.rstrip()
            columns = line.split("\t")
            cigar = parse_cigar(columns[5])
            if columns[0] in reads:
                read = reads[columns[0]]
            else:
                read = r.Read(columns[0], (len(columns[9]) + sum(cigar['H'])))
                reads[columns[0]] = read
            if int(columns[1]) & 4 or int(columns[4]) < opts.mapq:
                continue
            segment = s.Segment(segmentID, columns[0], columns[1], columns[2], columns[3], columns[4], len(columns[9]))
            segment.parseCigar(cigar)
            if float(segment.pid) < opts.pid:
                continue
            read.addSegment(segment)
            segments[segmentID] = segment
            segmentID += 1


def print_vcf_header():
    print(textwrap.dedent("""\
                ##fileformat=VCFv4.1
##fileDate=$date
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=BND,Description="Breakend">
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Distance between the two genomic positions ( END - POS )">
##INFO=<ID=RT,Number=3,Type=Integer,Description="Number of the different read types ( 2d , template , complement )">
##INFO=<ID=GAP,Number=1,Type=Integer,Description="Median number of bases between the two segments of the SV, in case of an insertion this is the size of the insertion">
##INFO=<ID=MAPQ,Number=2,Type=Integer,Description="Median mapping quality of the two segments of the structural variant">
##INFO=<ID=PID,Number=2,Type=Float,Description="Median percentage identity to the reference of the two segments of the structural variant">
##INFO=<ID=PLENGTH,Number=2,Type=Float,Description="Median segment length percentage of the two segments of the structural variant">
##INFO=<ID=RLENGTH,Number=1,Type=Integer,Description="Median length of the total reads">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=SQ,Number=1,Type=Float,Description="Phred-scaled probability that this site is variant (non-reference in this sample)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HR,Number=2,Type=Integer,Description="Number of hanging variant reads">
##FILTER=<ID=SVcluster,Description="There are more than $opt{cluster} SVs in a window of $opt{window} on both sides">
##FILTER=<ID=GAP,Description="The median gap size is larger than $opt{gap} for non insertions">
##FILTER=<ID=MapQual,Description="The median mapping quality is less than $opt{mapqf}">
##FILTER=<ID=PID,Description="The PID of one of the segments is less than $opt{pidf}">
##FILTER=<ID=CIPOS,Description="The CIPOS is greater or less than $opt{ci}">
##FILTER=<ID=CIEND,Description="The CIEND is greater or less than $opt{ci}">"""))
    print("\t".join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']))
    # print("threads", opts.threads)
    # print("sambamba", opts.sambamba)
    # print("split", opts.split)
    # print("pid", opts.pid)
    # print("mapq", opts.mapq)
    # print("distance", opts.distance)
    # print("count", opts.count)
    # print("refdistance", opts.refdistance)
    # print("unmapped", opts.unmapped)
    # print("matedistance", opts.matedistance)
    # print("window", opts.window)
    # print("cluster", opts.cluster)
    # print("mapqf", opts.mapqf)
    # print("pidf", opts.pidf)
    # print("gap", opts.gap)
    # print("ci", opts.ci)


# def check_opt():
#     print(sys.argv[1])
#     try:
#         if not re.match("/.bam$/", sys.argv[1]):
#             usage("No bam file given")
#             sys.exit()
#     except IndexError:
#         usage("No file given")
#         sys.exit()
#
#
# def usage(message):
#     print(message)


def main():
    # check_opt()
    print_vcf_header()
    parse_bam(bam, segmentID);
    parse_reads()
    parse_breakpoints()
    parse_svs()
    print(time.strftime("%c") + " Done")


if __name__ == '__main__':
    main()
