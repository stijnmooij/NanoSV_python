
class Variant:
    def __init__(self, chr, pos):
        self.chr = chr
        self.pos = pos
        self.segments = {}

    def add_segment(self, segment_id, variant_info):
        self.segments[segment_id[2]] = variant_info

    # def add_segment(self, segment_id, variant_info):
    #     self.segments[segment_id] = variant_info