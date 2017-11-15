# -*- coding: utf-8 -*-
import vcf
import re
import glob
import os.path
import hashlib

param = re.compile("LOW|MODIFIER")


class Vcf:
    def __init__(self, data):
        # 'CHROM', 'POS', 'ID', 'REF', 'ALT',
        # 'QUAL', 'FILTER', 'INFO', 'FORMAT',
        # 'start', 'end', 'alleles', 'samples',
        # '_sample_indexes', 'affected_start',
        # 'affected_end'
        # INFO:'AC', 'AF', 'AN', 'BaseQRankSum',
        # 'ClippingRankSum', 'DP', 'FS', 'MLEAC',
        # 'MLEAF', 'MQ', 'MQRankSum', 'QD',
        # 'ReadPosRankSum', 'SNPEFF_AMINO_ACID_CHANGE',
        # 'SNPEFF_CODON_CHANGE', 'SNPEFF_EFFECT',
        # 'SNPEFF_EXON_ID', 'SNPEFF_FUNCTIONAL_CLASS',
        # 'SNPEFF_GENE_BIOTYPE', 'SNPEFF_GENE_NAME',
        # 'SNPEFF_IMPACT', 'SNPEFF_TRANSCRIPT_ID', 'SOR'
        self.info = data.INFO
        self.alt = str(data.ALT[0])
        self.ref = data.REF
        self.pos = str(data.POS)
        self.annotation = data.INFO.get("SNPEFF_FUNCTIONAL_CLASS", "")
        self.impact = data.INFO.get("SNPEFF_IMPACT", "")
        self.gene = data.INFO.get("SNPEFF_GENE_NAME", "")
        self.feature = data.INFO.get("SNPEFF_TRANSCRIPT_ID", "")
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
    """read vcf file. only annotated"""
    filename = filepath(filename)
    result = []
    fp = open(filename, "r")
    vcf_r = vcf.Reader(fp)
    for v in vcf_r:
        vcfobj = Vcf(v)
        if vcfobj.impact is "":
            continue
        if re.match(param, vcfobj.impact) is None:
            result.append(vcfobj)
    fp.close()
    return result


def reference_read(filename):
    """read reference file"""
    filename = filepath(filename)
    fp = open(filename, "r")
    next(fp)
    reference = fp.read().strip().replace("\n", "")
    return reference


def provean_read(filename):
    f = open(filepath(filename), "r")
    gene = ""
    result = dict()
    p1 = re.compile("^# Query")
    p2 = re.compile("^#|^\[")
    for i in f:
        i = i.strip()
        if re.match(p1, i) is not None:
            gene = os.path.basename(i.split(":")[1].strip()).split(".")[0]
        if re.match(p2, i) is not None:
            continue
        scores = i.split()
        tmp = result.get(gene, {})
        tmp.update({scores[0]: scores[1]})
        result[gene] = tmp
    return result


def provean_union_add(provean_data, uniondata):
    result = dict()
    for gene, pscores in provean_data.items():
        vobjs = []
        for vobj in uniondata[gene]:
            vchange = vobj.info["SNPEFF_AMINO_ACID_CHANGE"]
            if (vchange in pscores) is False:
                continue
            vobj.provean_score = pscores[vchange]
            if float(vobj.provean_score) <= -2.5:
                vobj.provean_judge = "damaging"
            else:
                vobj.provean_judge = "tolerated"
            vobjs.append(vobj)
        if len(vobjs) > 0:
            result[gene] = vobjs
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


def gff_add(variants, gffdict=None):
    vsum = dict()
    for v in variants:
        if in_cobj(vsum.get(v.gene, []), v) is False:
            gff = gffdict[v.gene]
            v.start = gff.start
            v.end = gff.end
            v.strand = gff.strand
            v.description = gff.description
            vsum[v.gene] = vsum.get(v.gene, []) + [v]
    return vsum


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


def v_union(variant_list, gffdict=None):
    vsum = dict()
    if gffdict == None:
        for g in variant_list:
            for v in g:
                if in_cobj(vsum.get(v.gene, []), v) == False:
                    vsum[v.gene] = vsum.get(v.gene, []) + [v]
    else:
        for g in variant_list:
            for v in g:
                if in_cobj(vsum.get(v.gene, []), v) == False:
                    gff = gffdict[v.gene]
                    v.start = gff.start
                    v.end = gff.end
                    v.strand = gff.strand
                    v.description = gff.description
                    vsum[v.gene] = vsum.get(v.gene, []) + [v]
    return vsum


def in_cobj(listobj, cobj):
    tf = 0
    for i in listobj:
        if i.sha1 == cobj.sha1:
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
