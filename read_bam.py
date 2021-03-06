import os
import random
import sys
import textwrap
import time
import configparser

from collections import defaultdict

import read as r
import segment as s
import breakpoint as b
import sv as svclass

import re

cfg = configparser.ConfigParser()
cfg.read(sys.argv[2])

opts_threads = int(cfg.get('General Options', 'threads'))
opts_sambamba = cfg.get('General Options', 'sambamba')

opts_split = int(cfg.get('Filter Options', 'split'))
opts_pid = float(cfg.get('Filter Options', 'pid'))
opts_mapq = int(cfg.get('Filter Options', 'mapq'))

opts_distance = int(cfg.get('Detection Options', 'distance'))
opts_count = int(cfg.get('Detection Options', 'count'))
opts_refdistance = int(cfg.get('Detection Options', 'refdistance'))
opts_unmapped = int(cfg.get('Detection Options', 'unmapped'))
opts_matedistance = int(cfg.get('Detection Options', 'matedistance'))
opts_bed = cfg.get('Detection Options', 'bed')
opts_coverage_dupdel_check = cfg.get('Detection Options', 'coverage_dupdel_check')

opts_window = int(cfg.get('Output Filter Options', 'window'))
opts_cluster = int(cfg.get('Output Filter Options', 'cluster'))
opts_mapqf = int(cfg.get('Output Filter Options', 'mapqf'))
opts_pidf = float(cfg.get('Output Filter Options', 'pidf'))
opts_gap = int(cfg.get('Output Filter Options', 'gap'))
opts_ci = int(cfg.get('Output Filter Options', 'ci'))
opts_output = cfg.get('Output Filter Options', 'output')

bam = sys.argv[1]

reads = {}
segments = {}
breakpoints = {}
breakpoints_region = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
hanging_breakpoints_region = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
structural_variants = {}
structural_variants_region = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
structural_variants_region_2 = defaultdict(lambda: defaultdict(int))
output = open(opts_output, "w")
coverages = []
svID = 1
segmentID = 1


def parse_svs():
    global output
    sys.stderr.write(time.strftime("%c") + " Busy with reference reads...\n")
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
                    if sv_min > (segments[segment_id].pos + opts_refdistance) and \
                            sv_max < (segments[segment_id].end - opts_refdistance):
                        structural_variants[sv_id].format['DR'][x] += 1
                        structural_variants[sv_id].format['RO'][x] += (1 - 10 ** (-segments[segment_id].mapq / 10.0))
        prev_rname = segments[segment_id].rname

    sys.stderr.write(time.strftime("%c") + " Busy with printing to vcf...\n")
    for sv_id in sorted(structural_variants):
        if structural_variants[sv_id].set != 1: structural_variants[sv_id].setArguments()
        for sv_id_2 in sorted(structural_variants_region_2[structural_variants[sv_id].chr]):
            structural_variants[sv_id].setCluster(1)
            structural_variants[sv_id_2].setCluster(1)
            if sv_id_2 <= sv_id: continue
            if structural_variants[sv_id_2].set != 1: structural_variants[sv_id_2].setArguments()
            if structural_variants[sv_id].chr != structural_variants[sv_id_2].chr: continue
            if structural_variants[sv_id].chr2 != structural_variants[sv_id_2].chr2: continue
            if sv_id == 2370:
                print(abs(structural_variants[sv_id].pos - structural_variants[sv_id_2].pos) <= opts_window)
                print(abs(structural_variants[sv_id].info['END'] - structural_variants[sv_id_2].info['END']) <= opts_window)
            if abs(structural_variants[sv_id].pos - structural_variants[sv_id_2].pos) <= opts_window and \
                    abs(structural_variants[sv_id].info['END'] - structural_variants[sv_id_2].info['END']) <= opts_window:
                structural_variants[sv_id].SVcluster += 1
                structural_variants[sv_id_2].SVcluster += 1
            if structural_variants[sv_id_2].info['SVTYPE'] != "BND": continue
            if structural_variants[sv_id].info['SVTYPE'] != "BND": continue
            if re.match("\](\w+):(\d+)\]\w+", structural_variants[sv_id].alt) and not \
                    re.match("\w+\[(\w+):(\d+)\[", structural_variants[sv_id_2].alt): continue
            if re.match("\[(\w+):(\d+)\[\w+", structural_variants[sv_id].alt) and not \
                    re.match("\w+\](\w+):(\d+)\]", structural_variants[sv_id_2].alt): continue
            if re.match("\w+\](\w+):(\d+)\]", structural_variants[sv_id].alt) and not \
                    re.match("\[(\w+):(\d+)\[\w+", structural_variants[sv_id_2].alt): continue
            if re.match("\w+\[(\w+):(\d+)\[", structural_variants[sv_id].alt) and not \
                    re.match("\](\w+):(\d+)\]\w+", structural_variants[sv_id_2].alt): continue
            if abs(structural_variants[sv_id].pos - structural_variants[sv_id_2].pos) > opts_matedistance: continue
            if abs(structural_variants[sv_id].info['END'] - structural_variants[sv_id_2].info[
                'END']) > opts_matedistance: continue
            structural_variants[sv_id].info['MATEID'] = sv_id_2
            structural_variants[sv_id_2].info['MATEID'] = sv_id
        if sv_id == 2370:
            print(structural_variants[sv_id].SVcluster, opts_cluster, structural_variants[sv_id].filter)
        if structural_variants[sv_id].SVcluster > opts_cluster:
            structural_variants[sv_id].filter.append("SVcluster")
        if sv_id == 2370:
            print(structural_variants[sv_id].filter)
        if structural_variants[sv_id].info['GAP'] > opts_gap and structural_variants[sv_id].info['SVTYPE'] != "INS":
            structural_variants[sv_id].filter.append("GAP")
        if re.match("(\d+),(\d+)", structural_variants[sv_id].info['MAPQ']):
            ma = re.search("(\d+),(\d+)", structural_variants[sv_id].info['MAPQ'])
            if int(ma.group(1)) < opts_mapqf or int(ma.group(2)) < opts_mapqf:
                structural_variants[sv_id].filter.append("MapQual")
        if re.match("(\d.\d+),(\d.\d+)", structural_variants[sv_id].info['PID']):
            ma = re.search("(\d.\d+),(\d.\d+)", structural_variants[sv_id].info['PID'])
            if float(ma.group(1)) < opts_pidf or float(ma.group(2)) < opts_pidf:
                structural_variants[sv_id].filter.append("PID")
        if re.match("(\d+),(\d+)", structural_variants[sv_id].info['CIPOS']):
            ma = re.search("(\d+),(\d+)", structural_variants[sv_id].info['CIPOS'])
            if int(ma.group(1)) > opts_ci or int(ma.group(2)) > opts_ci:
                structural_variants[sv_id].filter.append("CIPOS")
        if re.match("(\d+),(\d+)", structural_variants[sv_id].info['CIEND']):
            ma = re.search("(\d+),(\d+)", structural_variants[sv_id].info['CIEND'])
            if int(ma.group(1)) > opts_ci or int(ma.group(2)) > opts_ci:
                structural_variants[sv_id].filter.append("CIEND")
        if structural_variants[sv_id].info['SVTYPE'] == 'INS':
            structural_variants[sv_id].info['SVLEN'] = structural_variants[sv_id].info['GAP']
        new_output = structural_variants[sv_id].printVCF(output)
        output = new_output


def addSVInfo(sv):
    if sum(sv.format['DV']) >= opts_count * 2:
        for flag in [sv.flag1 - 1, sv.flag2 + 1]:

            if flag == -1 or flag == 17: head_tail = 'T'
            if flag == 15 or flag == 1: head_tail = 'H'

            for hanging_pos in sorted(hanging_breakpoints_region[sv.chr][head_tail]):
                if hanging_pos < (min(sv.pos) - opts_distance):
                    continue
                if hanging_pos > (max(sv.pos) + opts_distance):
                    break
                for hanging_id in hanging_breakpoints_region[sv.chr][head_tail][hanging_pos]:
                    segment_id = hanging_breakpoints_region[sv.chr][head_tail][hanging_pos][hanging_id]
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
        structural_variants[sv.id] = sv
        structural_variants_region[sv.chr][min(sv.pos)][max(sv.pos)][sv.id] = 0
        structural_variants_region[sv.chr2][min(sv.info['END'])][max(sv.info['END'])][sv.id] = 1
        structural_variants_region_2[sv.chr][sv.id] = 1


def calculate_coverage_bed():
    global coverages
    os.system("cat " + opts_bed + "| sort > " + opts_bed[:-3] + "sorted.bed")
    original_coordinates = {}
    chromosomes = []
    total_positions = 0
    with open(opts_bed[:-3] + "sorted.bed", "r") as bed:
        for line in bed:
            line = line.rstrip().split("\t")
            line = [int(x) for x in line]
            if line[0] in original_coordinates:
                original_coordinates[line[0]].append([line[1], line[2]])
            else:
                original_coordinates[line[0]] = [[line[1], line[2]]]
            chromosomes.append(line[0])
            total_positions += (line[2] - line[1])
    file = open("bedfile.bed", "w")
    for i in range(1000000):
        random_position = random.randint(1, total_positions)
        position_string = ""
        for chromosome in chromosomes:
            find = False
            for coordinates in original_coordinates[chromosome]:
                if random_position <= (coordinates[1] - coordinates[0]):
                    position_string += str(chromosome) + "\t" + str(coordinates[0] + random_position - 1) + "\t" + \
                                       str(coordinates[0] + random_position)
                    find = True
                    break
                else:
                    random_position -= (coordinates[1] - coordinates[0])
            if find:
                break
        file.write(position_string + "\n")
    file.close()
    cmd = opts_sambamba + " depth base --min-coverage=0 " + bam + " -L bedfile.bed " + " | awk '{if (NR!=1) print $3}'"
    with os.popen(cmd) as coverageOutput:
        for coverage in coverageOutput:
            if coverage != ""and coverage != "\n":
                coverages.append(int(coverage))


def calculate_coverage_avg():
    global coverages
    with os.popen(opts_sambamba + ' view -s 0.10 ' + bam) as sample:
        for line in sample:
            line = line.rstrip()
            columns = line.split("\t")
            position = str(columns[2]) + ":" + str(columns[3]) + "-" + str(int(columns[3]) + 1)
            with os.popen(opts_sambamba + " depth base " + bam + " -L " + position + " | awk '{if (NR!=1) print $0}'") \
                    as coverageOutput:
                for coverage in coverageOutput:
                    coverages.append(coverage[2])


def parse_breakpoints_2(breakpoints_region_2):
    global coverages
    prev_pos_2 = -1
    global svID
    for pos_2 in sorted(breakpoints_region_2):
        for pos_1 in sorted(breakpoints_region_2[pos_2]):
            for breakpoint_id in breakpoints_region_2[pos_2][pos_1]:
                breakpoint = breakpoints[breakpoint_id]
                if (prev_pos_2 == -1):
                    sv = svclass.SV(svID, breakpoint, bam, opts_sambamba, coverages)
                    svID += 1
                elif abs(breakpoint.segment_2["pos"] - prev_pos_2) <= opts_distance:
                    sv.addBreakpoint(breakpoint)
                else:
                    addSVInfo(sv)
                    sv = svclass.SV(svID, breakpoint, bam, opts_sambamba, coverages)
                    svID += 1

                sv.addInfoField("PID",
                                [[segments[breakpoint.segment_1["id"]].pid],
                                 [segments[breakpoint.segment_2["id"]].pid]])
                sv.addInfoField("MAPQ",
                                [[segments[breakpoint.segment_1["id"]].mapq],
                                 [segments[breakpoint.segment_2["id"]].mapq]])
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
                prev_pos_2 = int(pos_2)
    addSVInfo(sv)


def parse_breakpoints():
    sys.stderr.write(time.strftime("%c") + " Busy with parsing breakpoints...\n")
    for region in breakpoints_region:
        prev_pos_1 = -1
        breakpoints_region_2 = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        for pos_1 in sorted(breakpoints_region[region]):
            for breakpoint_id in breakpoints_region[region][pos_1]:
                if prev_pos_1 == -1:
                    breakpoints_region_2[breakpoints[breakpoint_id].segment_2["pos"]][pos_1][breakpoint_id] = 1
                elif abs(pos_1 - prev_pos_1) <= opts_distance:
                    breakpoints_region_2[breakpoints[breakpoint_id].segment_2["pos"]][pos_1][breakpoint_id] = 1
                else:
                    parse_breakpoints_2(breakpoints_region_2)
                    breakpoints_region_2.clear()
                    breakpoints_region_2[breakpoints[breakpoint_id].segment_2["pos"]][pos_1][breakpoint_id] = 1
                prev_pos_1 = int(pos_1)
        parse_breakpoints_2(breakpoints_region_2)


def parse_reads():
    sys.stderr.write(time.strftime("%c") + " Busy with parsing read segments...\n")
    breakpointID = 1
    hanging_breakpointID = -1
    for qname in reads:
        clips = sorted(reads[qname].segments)
        if len(clips) == 1 or len(clips) > opts_split:
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
            if i == 0 and segment_1.clip >= opts_unmapped:
                hanging_breakpoint_pos = segment_1.pos
                if segment_1.flag == 16:
                    hanging_breakpoint_pos = segment_1.end
                    hanging_breakpoints_region[segment_1.rname]['T'][hanging_breakpoint_pos][
                        hanging_breakpointID] = segment_1.id
                    hanging_breakpointID -= 1
                else:
                    hanging_breakpoints_region[segment_1.rname]['H'][hanging_breakpoint_pos][
                        hanging_breakpointID] = segment_1.id
                    hanging_breakpointID = hanging_breakpointID - 1
            if i2 == len(clips) - 1 and segment_2.clip_2 >= opts_unmapped:
                hanging_breakpoint_pos = segment_2.end
                if segment_2.flag == 16:
                    hanging_breakpoint_pos = segment_2.pos
                    hanging_breakpoints_region[segment_2.rname]['H'][hanging_breakpoint_pos][
                        hanging_breakpointID] = segment_2.id
                    hanging_breakpointID -= 1
                else:
                    hanging_breakpoints_region[segment_2.rname]['T'][hanging_breakpoint_pos][
                        hanging_breakpointID] = segment_2.id
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
    sys.stderr.write(time.strftime("%c") + " Busy with parsing bam file...\n")

    with os.popen(opts_sambamba + ' view -t ' + str(opts_threads) + ' ' + bam) as bam:
        for line in bam:
            line = line.rstrip()
            columns = line.split("\t")
            cigar = parse_cigar(columns[5])
            if columns[0] in reads:
                read = reads[columns[0]]
            else:
                read = r.Read(columns[0], (len(columns[9]) + sum(cigar['H'])))
                reads[columns[0]] = read
            if int(columns[1]) & 4 or int(columns[4]) < opts_mapq:
                continue
            segment = s.Segment(segmentID, columns[0], columns[1], columns[2], columns[3], columns[4], len(columns[9]))
            segment.parseCigar(cigar)
            if float(segment.pid) < opts_pid:
                continue
            read.addSegment(segment)
            segments[segmentID] = segment
            segmentID += 1


def print_vcf_header():
    sample_name = os.system(opts_sambamba + ' view -H ' + bam + '| egrep "@RG" | egrep "SM:(\w+)"')
    if sample_name == "":
        sample_name = re.sub('(\.sorted)?\.bam$', '', str(bam))
    output.write(textwrap.dedent("""\
##fileformat=VCFv4.1
##fileDate=""" + str(time.strftime("%c")) + """
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
##INFO PVALue
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DR,Number=2,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=2,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=SQ,Number=1,Type=Float,Description="Phred-scaled probability that this site is variant (non-reference in this sample)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HR,Number=2,Type=Integer,Description="Number of hanging variant reads">
##FILTER=<ID=SVcluster,Description="There are more than """ + str(opts_cluster) + """ SVs in a window of """ + str(
        opts_window) + """ on both sides">
##FILTER=<ID=GAP,Description="The median gap size is larger than """ + str(opts_gap) + """ for non insertions">
##FILTER=<ID=MapQual,Description="The median mapping quality is less than """ + str(opts_mapqf) + """">
##FILTER=<ID=PID,Description="The PID of one of the segments is less than """ + str(opts_pidf) + """">
##FILTER=<ID=CIPOS,Description="The CIPOS is greater or less than """ + str(opts_ci) + """">
##FILTER=<ID=CIEND,Description="The CIEND is greater or less than """ + str(opts_ci) + """">\n"""))
    output.write("\t".join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', str(sample_name), "\n"]))



def main():
    print_vcf_header()
    if opts_coverage_dupdel_check == 'true':
        if opts_bed == '-':
            calculate_coverage_avg()
        else:
            calculate_coverage_bed()
    parse_bam(bam, segmentID);
    parse_reads()
    parse_breakpoints()
    parse_svs()
    output.close()
    sys.stderr.write(time.strftime("%c") + " Done\n")


if __name__ == '__main__':
    main()
