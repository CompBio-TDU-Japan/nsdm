#!/usr/bin/env python3

from . import fileparse
from . import seq
import re
import math


def allcheck(ref_file, vcfdir, gff_file):
    gff_dict = fileparse.gff_read(gff_file)
    reference = seq.Ref(ref_file)
    variants = fileparse.vcf_allread(vcfdir)
    summary = fileparse.v_union(variants, gff_dict)

    header = "\t".join([
        "gene",
        "ref",
        "alt",
        "pos_Pseq",
    ])
    print(header)
    for k, v in summary.items():
        reference.variant = v
        variantlist = reference.provean()
        if len(variantlist) > 0:
            for vinfo in variantlist:
                if vinfo.pref == vinfo.palt:
                    print(vinfo.gene,
                          "\t",
                          vinfo.pref,
                          "\t",
                          vinfo.palt,
                          "\t",
                          vinfo.pvp,
                          )
