
class Variant:
    def __init__(self, chr, pos):
        self.chr = chr
        self.pos = pos
        self.segments = {}

    def add_segment(self, segment_ID, variant_info):
        self.segments[segment_ID] = variant_info