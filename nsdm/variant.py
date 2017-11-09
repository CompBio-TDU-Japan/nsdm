#!/usr/bin/env python3

from . import fileparse


def allcheck(ref_file, vcfdir, gff_file):
    gff_dict = fileparse.gff_read(gff_file)
    variants = fileparse.vcf_allread(vcfdir)
    summary = fileparse.v_union(variants, gff_dict)

    header = "\t".join([
        "gene",
        "aa_change",
    ])
    print(header)
    for k, variantlist in summary.items():
        for v in variantlist:
            if v.info["SNPEFF_FUNCTIONAL_CLASS"] != "MISSENSE":
                continue
            aach = v.info["SNPEFF_AMINO_ACID_CHANGE"]
            payload = "\t".join([
                v.gene,
                aach,
            ])
            if aach[0] == aach[-1]:
                print("#same!" + payload)
            else:
                print(payload)
                print(payload, v.info["SNPEFF_EFFECT"])
#    for k, v in summary.items():
#        reference.variant = v
#        variantlist = reference.provean()
#        if len(variantlist) > 0:
#            for vinfo in variantlist:
#                if vinfo.pref == vinfo.palt:
#                    print("\t".join([
#                        vinfo.gene,
#                        str(vinfo.pvp),
#                        vinfo.pref,
#                        vinfo.palt,
#                        str(vinfo.pos),
#                        str(vinfo.ref),
#                        str(vinfo.alt),
#                    ]))
