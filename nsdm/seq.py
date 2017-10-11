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
            if v.annotation != "missense_variant":
                continue
            pos = int(v.pos)
            vnseq[pos - 1] = v.alt
            v.nvp = pos - (int(start) + 1)
            v.pvp = math.ceil((v.nvp + 1) / 3) - 1
            if v.strand == "-":
                v.nvp = len(seq) - (v.nvp) - 1
                v.pvp = math.ceil((v.nvp + 1) / 3) - 1
                print(f"{v.strand}|", v.gene, v.pos)
                print(v.ref, "=>", v.alt)
                print(seq_reverse(seq)[v.nvp], "=>",
                      seq_reverse(vnseq[start:end])[v.pvp])
            print(f"{v.strand}|", v.gene, v.pos)
            print(v.ref, "=>", v.alt)
            print(seq[v.nvp], "=>", vnseq[pos - 1])
            exit()
            result.append(v)
        vseq = "".join(vnseq)[start:end]
        if self.variant[0].strand == "-":
            seq = seq_reverse(seq)
            vseq = seq_reverse(vseq)
        nseq = seq
        nvseq = vseq
        vppos = [x.pvp for x in result]
        vinfov = []
        vinfon = []
        if self.variant[0].strand == "-":
            nseq = seq_reverse(seq)
            nvseq = seq_reverse(vseq)
            seq, vinfon = translate(nseq, vppos)
            vseq, vinfov = translate(nvseq, vppos)
        else:
            seq, vinfon = translate(seq, vppos)
            vseq, vinfov = translate(vseq, vppos)
        if result[0].strand == "-":
            for n, v in enumerate(result):
                print("--")
                print(len(seq))
                print(v.ref, "=>", v.alt, v.pvp)
                print(len(seq), "/", seq[v.pvp], "=>",
                      len(seq), "/", vseq[v.pvp])
                v.palt = vseq[v.pvp]
                v.pref = seq[v.pvp]
                v.nseq = nseq
                v.nvseq = nvseq
                v.change = [x for x in zip(vinfon, vinfov)]
                print(v.change)
                result[n] = v
            exit()
        return (seq.split("*")[0], result)


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
                    variants.append(s + "|" + str(n))
                ret = ret + p
                break
            elif re.match(pattern, s):
                if n in variant:
                    variants.append(s + "|" + str(n))
                ret = ret + "X"
                break
    return ret, variants
