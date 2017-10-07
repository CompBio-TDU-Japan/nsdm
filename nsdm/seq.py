#!/usr/bin/env python3

from . import fileparse


class Ref:
    def __init__(self, reference_file):
        self.seq = fileparse.reference_read(reference_file)

    def cut(self, start, end):
        return self.seq[start - 1:end]


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
