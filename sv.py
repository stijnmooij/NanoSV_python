#!/usr/bin/python

from statistics import median
from math import log10
from math import floor


import math
import re

import os


class SV:
    def __init__(self, id, breakpoint, bam, sambamba, avg_sample_coverage):
        self.bam = bam
        self.sambamba = sambamba
        self.avg_sample_coverage = avg_sample_coverage
        self.id = id
        self.chr = breakpoint.segment_1["rname"]
        self.chr2 = breakpoint.segment_2["rname"]
        self.flag1 = breakpoint.segment_1["flag"]
        self.flag2 = breakpoint.segment_2["flag"]
        self.pos = [breakpoint.segment_1["pos"]]
        self.ref = 'N'
        self.alt = "<" + str(breakpoint.svtype) + ">"
        self.qual = '.'
        self.filter = ['PASS']
        self.info = {
            'PRECISE': "IMPRECISE",
            'END': [breakpoint.segment_2["pos"]],
            'SVTYPE': breakpoint.svtype,
            'SVMETHOD': 'nanosv',
        }
        self.format = {
            'GT': './.',
            'DV': [1, 1],
            'VO': [(1 - 10 ** (-breakpoint.segment_1["mapq"] / 10.0)), (1 - 10 ** (-breakpoint.segment_2["mapq"] / 10.0))],
            'DR': [0, 0],
            'RO': [0, 0],
            'HR': [0, 0]
        }
        self.breakpoints = [breakpoint.id]
        self.set = 0

    def addBreakpoint(self, breakpoint):
        self.breakpoints.append(breakpoint.id)
        self.pos.append(breakpoint.segment_1["pos"])
        self.info['END'].append(breakpoint.segment_2["pos"])
        self.format['DV'][0] += 1
        self.format['DV'][1] += 1
        self.format['VO'][0] += (1 - 10 ** (-breakpoint.segment_1["mapq"] / 10.0))
        self.format['VO'][1] += (1 - 10 ** (-breakpoint.segment_2["mapq"] / 10.0))

    def addInfoField(self, key, value):
        if key in self.info:
            if isinstance(self.info[key], list):
                try:
                    if isinstance(self.info[key][0], list) or isinstance(self.info[key][1], list):
                        if value[0] is not None: self.info[key][0].append(value[0][0])
                        if value[1] is not None: self.info[key][1].append(value[1][0])
                    else:
                        self.info[key].append(value[0])

                except IndexError:
                    self.info[key].append(value[0])
        else:
            self.info[key] = value

    def significanceTest(self, avg_dupdel_cov, dup):
        teller = 0
        if dup:
            for value in self.avg_sample_coverage:
                if float(value) > avg_dupdel_cov:
                    teller += 1
        else:
            for value in self.avg_sample_coverage:
                if float(value) < avg_dupdel_cov:
                    teller += 1
        if teller/len(self.avg_sample_coverage) < 0.05:
            return True
        else:
            return False

    def setArguments(self):
        self.info['CIPOS'] = str(int(min(self.pos) - median(self.pos))) + "," + str(int(max(self.pos) - median(self.pos)))
        self.info['CIEND'] = str(int(min(self.info['END']) - median(self.info['END']))) + "," + str(int(
            max(self.info['END']) - median(self.info['END'])))
        if self.info['CIPOS'] == "0,0" and self.info['CIEND']: self.info['PRECISE'] = "PRECISE"
        self.pos = median(self.pos)
        self.info['END'] = floor(median(self.info['END']))
        self.info['SVLEN'] = (self.info['END'] - self.pos)
        self.setInfoField()

        dup = 0
        if self.alt == "<BND>":
            if self.flag1 == 16:
                if self.flag2 == 16:
                    self.alt = "]" + self.chr2 + ":" + str(floor(self.info['END'])) + "]" + self.ref
                    avg_dupdel_cov = self.getDupDelcoverage()
                    if self.significanceTest(avg_dupdel_cov, True):
                        self.alt = "<DUP>"
                    dup = 1
                else:
                    self.alt = "[" + self.chr2 + ":" + str(floor(self.info['END'])) + "[" + self.ref
            else:
                if self.flag2 == 16:
                    self.alt = self.ref + "]" + self.chr2 + ":" + str(floor(self.info['END'])) + "]"
                else:
                    self.alt = self.ref + "[" + self.chr2 + ":" + str(floor(self.info['END'])) + "["
                    avg_dupdel_cov = self.getDupDelcoverage()
                    if self.significanceTest(avg_dupdel_cov, False):
                        self.alt = "<DEL>"
        gt_lplist = self.bayes_gt(sum(self.format['RO']), sum(self.format['VO']), dup)
        gt_idx = gt_lplist.index(max(gt_lplist))

        gt_sum = 0
        for gt in gt_lplist:
            gt_sum += 10 ** gt

        if (gt_sum > 0):
            gt_sum_log = log10(gt_sum)
            sample_qual = abs(-10 * (gt_lplist[0] - gt_sum_log))

            if 1 - (10 ** gt_lplist[gt_idx] / 10 ** gt_sum_log) == 0:
                phred_gq = 200
            else:
                phred_gq = abs(-10 * log10(1 - (10 ** gt_lplist[gt_idx] / 10 ** gt_sum_log)))


            self.format['GQ'] = int(phred_gq)
            self.format['SQ'] = format(sample_qual, '.3f')

            if gt_idx == 0:
                self.format['GT'] = '0/0'
            elif gt_idx == 1:
                self.format['GT'] = '0/1'
            elif gt_idx == 2:
                self.format['GT'] = '1/1'

        self.set = 1

    def getDupDelcoverage(self):
        start = round(self.pos)
        stop = round(self.info['END'])
        if start > stop:
            start = round(self.info['END'])
            stop = round(self.pos)
        position = str(self.chr) + ":" + str(start) + "-" + str(stop)
        dupdel_coverages = []
        with os.popen(self.sambamba + " depth base " + self.bam + " -L " + position + " | awk '{if (NR!=1) print $3 }'") as commandoutput:
            for line in commandoutput:
                dupdel_coverages.append(float(line.rstrip()))
        return float(sum(dupdel_coverages))/len(dupdel_coverages)

    def log_choose(self, n, k):
        r = 0.0
        # swap for efficiency if k is more than half of n
        if k * 2 > n:
            k = n - k
        for d in range(1, k + 1):
            r += math.log(n, 10)
            r -= math.log(d, 10)
            n -= 1

        return (r)

    def bayes_gt(self, ref, alt, is_dup):
        # probability of seeing an alt read with true genotype of of hom_ref, het, hom_alt respectively
        if is_dup:  # specialized logic to handle non-destructive events such as duplications
            p_alt = [1e-2, 1/3.0, 0.5]
        else:
            p_alt = [1e-3, 0.5, 0.9]

        total = ref + alt
        log_combo = self.log_choose(int(total), int(alt))

        lp_homref = log_combo + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
        lp_het = log_combo + (alt * math.log(p_alt[1], 10)) + (ref * math.log(1 - p_alt[1], 10))
        lp_homalt = log_combo + alt * math.log(p_alt[2], 10) + ref * math.log(1 - p_alt[2], 10)

        return [lp_homref, lp_het, lp_homalt]

    def setInfoField(self):
        for field in self.info:
            if field == "RT":
                rt = [0, 0, 0]
                for type in self.info[field]:
                    if type == "2d": rt[0] += 1
                    if type == "template": rt[1] += 1
                    if type == "complement": rt[2] += 1
                self.info[field] = ",".join(map(str, rt))
            elif isinstance(self.info[field], list):
                if isinstance(self.info[field][0], list):
                    self.info[field][0] = median(map(float, self.info[field][0]))
                    self.info[field][1] = median(map(float, self.info[field][1]))
                    if field == "MAPQ":
                        self.info[field][0] = int(self.info[field][0])
                        self.info[field][1] = int(self.info[field][1])
                    elif field == "PID" or field == "PLENGTH":
                        self.info[field][0] = round(self.info[field][0], 3)
                        self.info[field][1] = round(self.info[field][1], 3)
                    self.info[field] = ",".join(map(str, self.info[field]))
                else:
                    self.info[field] = self.median_type(self.info[field])

    def setCluster(self, SVcluster):
        self.SVcluster = SVcluster

    def printVCF(self, output):
        if len(self.filter) == 1:
            self.filter = "PASS"
        else:
            self.filter = ";".join(self.filter)
            self.filter = self.filter.replace('PASS;', '')

        output.write("\t".join(map(str, [self.chr, int(self.pos), self.id, self.ref, self.alt, self.qual, self.filter, self.info['PRECISE']])))
        for field in self.info:
            if field == 'PRECISE':
                continue
            if field == "END" and self.info['SVTYPE'] == "BND":
                continue
            if field == "SVLEN" and not re.match("/^$self->{_chr2}$/", self.chr):
                continue

            output.write(";%s=%s" % (field, self.info[field]))
        output.write("\tGT")
        values = []
        for field in self.format:
            if field == 'GT':
                continue
            if not re.match("DR|DV|HR|SQ|GQ", field):
                continue
            output.write(":%s" % (field))
            value = self.format[field]
            if isinstance(value, list):
                value = ",".join(map(str, self.format[field]))
            values.append(str(value))
        output.write("\t%s:%s\n" % (self.format['GT'], ":".join(values)))
        return output

    def median_type(self, toCalculate):
        if type(toCalculate[0]) == float:
            return median(toCalculate)
        elif type(toCalculate[0]) == int:
            return int(median(toCalculate))
