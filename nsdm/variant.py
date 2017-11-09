#!/usr/bin/env python3

from . import fileparse


def allcheck(ref_file, vcfdir, gff_file):
    gff_dict = fileparse.gff_read(gff_file)
    variants = fileparse.vcf_allread(vcfdir)
    summary = fileparse.v_union(variants, gff_dict)

    result = dict()
    print("#", "\t".join(["gene", "change", "description"]))
    for k, variantlist in summary.items():
        description = variantlist[0].description
        for v in variantlist:
            if v.info["SNPEFF_EFFECT"] == "SYNONYMOUS_CODING":
                continue
            aach = v.info["SNPEFF_AMINO_ACID_CHANGE"]
            payload = "|".join([
                aach,
                v.annotation,
            ])
            result[k] = result.get(k, []) + [payload]
        print(k + "\t" + str(result[k]) + "\t" + description)
