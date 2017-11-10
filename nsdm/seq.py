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
        for variant in variantlist:
            print(variant.__dict__)
#        x = self.variant[0]
#        start = 0
#        end = 0
#        if isinstance(x.start, str):
#            start = int(x.start) - 1
#        if isinstance(x.end, str):
#            end = int(x.end)
#        genome = self.seq
#        seq = self.seq[start:end]
#        base_vpseq_genome = genome
#        base_vpseq_genome = list(base_vpseq_genome)
#        result = []
#        for v in self.variant:
#            if v.annotation != "missense_variant":
#                continue
#            pos = int(v.pos) - 1
#            base_vpseq_genome[pos] = v.alt
#            v.nvp = pos - start
#            v.pvp = math.ceil((v.nvp + 1) / 3) - 1
#            if v.strand == "-":
#                v.nvp = len(seq) - (v.nvp) - 1
#                v.pvp = math.ceil((v.nvp + 1) / 3) - 1
#            result.append(v)
#        vseq = "".join(base_vpseq_genome)[start:end]
#        if len(result) == 0:
#            return result
#        if result[0].strand == "-":
#            seq = seq_reverse(seq)
#            vseq = seq_reverse(vseq)
#        vppos = [x.pvp for x in result]
#        nseq = seq
#        nvseq = vseq
#        vinfov = []
#        vinfon = []
#        if result[0].strand == "-":
#            seq, vinfon = translate(seq, vppos)
#            vseq, vinfov = translate(vseq, vppos)
#        else:
#            seq, vinfon = translate(seq, vppos)
#            vseq, vinfov = translate(vseq, vppos)
#        for n, v in enumerate(result):
#            v.palt = vseq[v.pvp]
#            v.pref = seq[v.pvp]
#            v.nseq = nseq
#            v.nvseq = nvseq
#            v.protein = seq.split("*")[0]
#            v.vprotein = vseq.split("*")[0]
#            v.codon_aa = [x for x in zip(vinfon, vinfov)]
#            result[n] = v
#        return result
#


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
