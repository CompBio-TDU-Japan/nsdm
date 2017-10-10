#!/usr/bin/env python3

from . import fileparse
import re
import math


class Ref:
    def __init__(self, reference_file):
        self.seq = fileparse.reference_read(reference_file)
        self.variant = ""

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

    def provean(self):
        x = self.variant[0]
        start = 0
        end = 0
        if isinstance(x.start, str):
            start = int(x.start) - 1
        if isinstance(x.end, str):
            end = int(x.end)
        genome = self.seq
        seq = self.seq[start:end]
        vnseq = genome
        vnseq = list(vnseq)
        result = []
        for v in self.variant:
            pos = int(v.pos)
            vnseq[pos - 1] = v.alt
            v.nvp = pos - (int(v.start) + 1)
            v.pvp = math.ceil(v.nvp / 3) - 1
            if v.strand == "-":
                v.nvp = len(seq) - (v.nvp + 1)
                v.pvp = math.ceil(v.nvp / 3) - 1
            result.append(v)
        vseq = "".join(vnseq)[start:end]

        if self.variant[0].strand == "-":
            seq = translate(seq_reverse(seq))
            vseq = translate(seq_reverse(vseq))
        else:
            seq = translate(seq)
            vseq = translate(vseq)
        for n, v in enumerate(result):
            v.palt = vseq[v.pvp]
            v.pref = seq[v.pvp]
            result[n] = v
        return (seq.split("*")[0], result)


def seq_reverse(seq):
    compliments = {'N': 'N', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    ret = "".join([compliments[x] for x in seq])[::-1]
    return ret


def translate(seq):
    pattern = re.compile(r"N")
    AAs = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    Base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    Base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    Base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
    target = re.findall('.' * 3, seq)
    ret = ""
    for s in target:
        for (i1, i2, i3, p) in zip(Base1, Base2, Base3, AAs):
            if s == i1 + i2 + i3:
                ret = ret + p
                break
            elif re.match(pattern, s):
                ret = ret + "X"
                break
    return ret
