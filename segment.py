#!/usr/bin/python

class Segment:
    def __init__(self, id, qname, flag, rname, pos, mapq, length):
        self.id = id
        self.qname = qname
        self.flag = int(flag)
        self.rname = rname
        self.pos = int(pos)
        self.mapq = int(mapq)
        self.length = int(length)
        self.end = (int(pos) + int(length))
        self.clip = False
        self.clip_2 = False
        self.pid = False

    def parseCigar(self, cigar):
        if self.flag & 16:
            self.clip = cigar['H'][1]
            self.clip_2 = cigar['H'][0]
        else:
            self.clip = cigar['H'][0]
            self.clip_2 = cigar['H'][1]

        self.end += cigar['D']
        self.end -= cigar['I']
        self.pid = format(cigar['='] / float(self.length), '.3f')

    def setPlength(self, rlength):
        self.plength = format(self.length/ rlength, '.3f')