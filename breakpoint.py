#!/usr/bin/python

class Breakpoint:
    def __init__(self, id, segment_1, segment_2):
        self.id = id
        self.segment_1 = segment_1
        self.segment_2 = segment_2

    def setGap(self, gap):
        self.gap = gap

    def setBreakpoint(self, segment_1, segment_2):
        self.segment_1.pos = segment_1.pos if segment_1.flag == 16 else segment_1.end
        self.segment_2.pos = segment_2.pos if segment_2.flag == 16 else segment_2.end

    def switchSegments(self):
        self.segment_1, self.segment_2 = self.segment_2, self.segment_1
        self.segment_1.flag = 0 if self.segment_1.flag == 16 else 16
        self.segment_2.flag = 0 if self.segment_2.flag == 16 else 16

    def setSVtype(self):
        if self.segment_1.rname == self.segment_2.rname:
            self.svtype = "BND"
        else:
            if self.segment_1.flag == 16:
                if self.segment_2.flag == 16:
                    self.svtype = "DUP"
                else:
                    self.svtype = "BND"
            elif self.segment_2.flag == 16:
                self.svtype = "BND"
            else:
                self.svtype = "DEL"
            if abs(self.segment_2.pos - self.segment_1.pos) < self.gap:
                self.svtype = "INS"
                self.segment_2.pos = self.segment_1.pos+1
                self.segment_1.flag = 0
                self.segment_2.flag = 0
