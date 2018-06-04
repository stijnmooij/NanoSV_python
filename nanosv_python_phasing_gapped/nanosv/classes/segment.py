#!/usr/bin/python


class Segment:
    # def __init__(self, id, qname, flag, rname, pos, mapq, length, variations):
    #     self.id = id
    #     self.qname = qname
    #     self.flag = int(flag)
    #     self.rname = rname
    #     self.pos = int(pos)
    #     self.mapq = int(mapq)
    #     self.length = int(length)
    #     self.end = (int(pos) + int(length))
    #     self.clip = False
    #     self.clip_2 = False
    #     self.pid = False
    #     self.variations = variations
    def __init__(self, qname, flag, rname, pos, mapq, length, clip, clip_2, pid):
        #self.qname_clip = str(rname) + ":" + str(pos) + "-" + str(qname)
        self.id = [rname, pos, qname + ";" + str(clip)]
        self.qname = qname
        self.flag = int(flag)
        self.rname = rname
        self.pos = int(pos)
        self.mapq = int(mapq)
        self.length = int(length)
        self.end = (int(pos) + int(length))
        self.clip = clip
        self.clip_2 = clip_2
        self.pid = pid

    def setPlength(self, rlength):
        """
        Calculates plength
        :param rlength is used to calculate the plength:
        """
        self.plength = format(self.length / rlength, '.3f')
