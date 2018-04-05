import NanoSV
from utils import parse_breakpoints as breakpoint
from utils import create_vcf as vcf


def report_svs_to_vcf():
    for sv_id in sorted(breakpoint.structural_variants):
        breakpoint.structural_variants[sv_id].printVCF(vcf.vcf_writer)
