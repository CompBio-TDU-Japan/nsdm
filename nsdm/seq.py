#!/usr/bin/env python3

from . import fileparse
import re
import math


class Ref:
    def __init__(self, reference_file):
        self.seq = fileparse.reference_read(reference_file)

    def cut(self):
        x = self.variant[0]
        start = 0
        end = 0
        if isinstance(x.start, str):
            start = int(x.start) - 1
        if isinstance(x.end, str):
            end = int(x.end)
        seq = self.seq[start:end]
        vseq = self.seq
        vseq = list(vseq)
        for v in self.variant:
            pos = int(v.pos)
            vseq[(pos - 1)] = v.alt
        vseq = "".join(vseq)[start:end]
        if self.variant[0].strand == "-":
            seq = translate(seq_reverse(seq))
            vseq = translate(seq_reverse(vseq))
        else:
            seq = translate(seq)
            vseq = translate(vseq)
        return (seq.split("*")[0], vseq.split("*")[0])

    def provean(self, variantlist):
        x = variantlist[0]
        result = dict()
        gene = x.gene
        start = int(x.start) - 1
        end = int(x.end)
        var = []
        for v in variantlist:
            if v.annotation != "MISSENSE":
                continue
            var.append(v.info["SNPEFF_AMINO_ACID_CHANGE"])
        result["variant"] = var
        nseq = self.seq[start:end]
        if variantlist[0].strand == "-":
            nseq = seq_reverse(nseq)
        pseq, none = translate(nseq)
        result["fasta"] = [gene,
                           ">" + gene + "\n" + pseq.split("*")[0]]
        return result


def seq_reverse(seq):
    compliments = {'N': 'N', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    ret = "".join([compliments[x] for x in seq])[::-1]
    return ret


def translate(seq, variant=[]):
    pattern = re.compile(r"N")
    AAs = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    Base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    Base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    Base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
    target = re.findall('.' * 3, seq)
    ret = ""
    variants = []
    for n, s in enumerate(target):
        for (i1, i2, i3, p) in zip(Base1, Base2, Base3, AAs):
            if s == i1 + i2 + i3:
                if n in variant:
                    variants.append(s + "|" + str(n) + "|" + p)
                ret = ret + p
                break
            elif re.match(pattern, s):
                if n in variant:
                    variants.append(s + "|" + str(n))
                ret = ret + "X"
                break
    return (ret, variants)
