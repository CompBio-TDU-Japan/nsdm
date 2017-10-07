# -*- coding: utf-8 -*-
import vcf
import re
import glob
import os.path
import hashlib

param = re.compile("LOW|MODIFIER")


class Vcf:
    def __init__(self, data):
        info = data.INFO["ANN"][0].split("|")
        self.alt = str(data.ALT[0])
        self.ref = data.REF
        self.pos = str(data.POS)
        self.annotation = info[1]
        self.impact = info[2]
        self.gene = info[4]
        self.feature = info[6].split(".")[0]
        vid = (self.alt + self.ref + self.pos + self.annotation
               + self.impact + self.gene + self.feature)
        self.sha1 = hashlib.sha1(vid.encode('utf-8')).hexdigest()


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


def reference_read(filename):
    """read reference file"""
    filename = filepath(filename)
    fp = open(filename, "r")
    reference = "".join([x.strip() for x in fp.readlines()[1:]])
    return reference


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


def v_intersect(variant_list):
    vsum = dict()
    for g in variant_list:
        for v in g:
            vsum[v.sha1] = vsum.get(v.sha1, []) + [v]
    ret = []
    for k, v in vsum.items():
        if len(v) == len(variant_list):
            ret.append(v[0])
    return ret


def v_union(variant_list):
    vsum = dict()
    tf = 0
    for g in variant_list:
        for v in g:
            if in_cobj(vsum.get(v.gene, []), v) == False:
                vsum[v.gene] = vsum.get(v.gene, []) + [v]
    return vsum


def in_cobj(listobj, cobj):
    tf = 0
    for i in listobj:
        if i.__dict__ == cobj.__dict__:
            tf += 1
    if tf > 0:
        return True
    else:
        return False


def bbhdict(f):
    filename = filepath(f)
    f = open(filename, "r")
    ret = dict()
    for i in f:
        k, v = i.strip().split("\t")[:2]
        ret[k] = v
    return ret


def bbh(f, variantdict):
    result = dict()
    bdict = bbhdict(f)
    for key, valuelist in variantdict.items():
        homolog = bdict[key]
        if homolog == "-":
            continue
        for v in valuelist:
            v.homolog = homolog
            result[homolog] = result.get(homolog, []) + [v]
    return result


def filepath(file):
    return os.path.abspath(os.path.expanduser(file))
