# -*- coding: utf-8 -*-
import vcf


##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">

class Vcf:
    def __init__(self,data):
        info = data.INFO["ANN"][0].split("|")
        self.alt = data.REF
        self.ref = data.REF
        self.pos = data.POS
        self.annotation = info[1]
        self.impact = info[2]
        self.gene = info[4]
        self.feature = info[6].split(".")[0]

class Gff:
    def __init__(self,data):
        self.start = data[3]
        self.end = data[4]
        self.strand = data[6]
        info = data[8].split(";")
        self.gene = info[0].split(":")[1]
        self.description = info[2].split("=")[1]


def vcf_read(filename):
    """read vcf file"""
    result = []
    fp = open(filename,"r")
    vcf_r = vcf.Reader(fp)
    for v in vcf_r:
        result.append(Vcf(v))
    fp.close()
    return result

def gff_read(filename):
    """read gff file"""
    result = dict()
    fp = open(filename,"r")
    for i in fp:
        data = i.split()
        if i[0] != "#" and data[2] == "gene":
            gene = data[8].split(";")[0].split(":")[1]
            result[gene] = Gff(data)
    fp.close()
    return result
