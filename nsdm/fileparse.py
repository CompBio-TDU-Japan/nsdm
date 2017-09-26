# -*- coding: utf-8 -*-
import vcf
import re
import glob
import os.path

param = re.compile("LOW|MODIFIER")


class Vcf:
    def __init__(self, data):
        info = data.INFO["ANN"][0].split("|")
        self.alt = data.REF
        self.ref = data.REF
        self.pos = data.POS
        self.annotation = info[1]
        self.impact = info[2]
        self.gene = info[4]
        self.feature = info[6].split(".")[0]


class Gff:
    def __init__(self, data):
        self.start = data[3]
        self.end = data[4]
        self.strand = data[6]
        info = data[8].split(";")
        self.gene = info[0].split(":")[1]
        self.description = info[2].split("=")[1].replace("%2C", ",")


def vcf_read(filename):
    """read vcf file"""
    filename = filepath(filename)
    result = []
    fp = open(filename, "r")
    vcf_r = vcf.Reader(fp)
    for v in vcf_r:
        vcfobj = Vcf(v)
        if re.match(param, vcfobj.impact) == None:
            result.append(vcfobj)
    fp.close()
    return result


def vcf_allread(targetdir):
    result = []
    allvariantdir = filepath(targetdir)
    vcffiles = sorted(glob.glob(allvariantdir + "/*.vcf"))
    for i in vcffiles:
        vlist = vcf_read(i)
        result.append(vlist)
    return result


def gff_read(filename):
    """read gff file"""
    filename = filepath(filename)
    result = dict()
    fp = open(filename, "r")
    for i in fp:
        data = i.split("\t")
        if i[0] != "#" and data[2] == "gene":
            gene = data[8].split(";")[0].split(":")[1]
            result[gene] = Gff(data)
    fp.close()
    return result

def filepath(file):
    return os.path.abspath(os.path.expanduser(file))
