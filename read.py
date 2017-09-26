#!/usr/bin/python

class Read:
    def __init__(self, qname, length):
        self.qname = qname
        self.length = length
        self.segments = dict()

    def addSegment(self, segment):
        self.segments[segment.clip] = segment.id


    # @property
    # def qname(self):
    #     """I'm the 'qname'property"""
    #     print("getter of qname called")
    #     return self.qname
    #
    # @property
    # def length(self):
    #     """I'm the 'length' property"""
    #     print("getter of length called")
    #     return self.length




# self.segments heeft segment als key en id als value


#vragen: wat doet clip precies?