# -*- coding: utf-8 -*-
import vcf
import re
import glob
import os.path
from Bio import SeqUtils as SeqUtils
import hashlib

param = re.compile("LOW|MODIFIER")


class Vcf:
    def __init__(self, data):
        # VCF object. The following are the attributes.
        # 'Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID',
        # 'Feature_Type', 'Feature_ID', 'Transcript_BioType',
        # 'Rank', 'HGVSc', 'HGVSp', 'cDNApos/cDNAlength',
        # 'CDSpos/CDSlength', 'AApos/AAlength', 'Distance',
        # 'ERRORS/WARNINGS/INFO', 'AC', 'AF', 'AN', 'BaseQRankSum',
        # 'ClippingRankSum', 'DP', 'FS',
        # 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR',
        # 'CHROM', 'POS', 'ID',
        # 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT', 'start', 'end', 'alleles',
        # 'samples', '_sample_indexes','sample_name'
        # 'affected_start', 'affected_end', 'sha1'
        annheader = "".join(['Allele|Annotation|Annotation_Impact|',
                             'Gene_Name|Gene_ID|Feature_Type|',
                             'Feature_ID|Transcript_BioType|Rank|',
                             'HGVSc|HGVSp|cDNApos/cDNAlength',
                             '|CDSpos/CDSlength|AApos/AAlength|',
                             'Distance|ERRORS/WARNINGS/INFO']).split("|")
        ann = {k: v for k, v in zip(
            annheader, data.INFO.get("ANN")[0].split("|"))}
        data.INFO.pop("ANN")
        self.__dict__.update(ann)
        self.__dict__.update(data.INFO)
        data.__dict__.pop("INFO")
        self.__dict__.update(data.__dict__)
        self.ALT = str(self.ALT[0])
        checkid = "".join([
            self.REF,
            self.ALT,
            self.Annotation_Impact,
            str(self.POS),
            self.Gene_ID,
            self.Feature_ID,
            self.Annotation,
        ])
        self.sha1 = hashlib.sha1(checkid.encode('utf-8')).hexdigest()
        self.sample_name = list(self._sample_indexes)[0]
        self.SNPEFF_AMINO_ACID_CHANGE = ""
        if self.Annotation == "missense_variant":
            aach = self.HGVSp.split(".")[1]
            aminochange = [aach[:3], aach[3:-3], aach[-3:]]
            self.SNPEFF_AMINO_ACID_CHANGE = "".join([
                SeqUtils.IUPACData.protein_letters_3to1_extended[aminochange[0]],
                aminochange[1],
                SeqUtils.IUPACData.protein_letters_3to1_extended[aminochange[2]],
            ])
        elif self.Annotation == "stop_gained":
            aach = self.HGVSp.split(".")[1]
            aminochange = [aach[:3], aach[3:-3], aach[-3:]]
            self.SNPEFF_AMINO_ACID_CHANGE = "".join([
                SeqUtils.IUPACData.protein_letters_3to1_extended[aminochange[0]],
                aminochange[1],
                aminochange[2].strip(),
            ])


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
        if vcfobj.Annotation_Impact is "":
            continue
        if re.match(param, vcfobj.Annotation_Impact) is None:
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


def provean_add(provean_data, datadict):
    result = dict()
    for gene, pscores in provean_data.items():
        vobjs = []
        if (gene in datadict) is False:
            continue
        for vobj in datadict[gene]:
            vchange = vobj.SNPEFF_AMINO_ACID_CHANGE
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


def provean_union_add(provean_data, uniondata):
    result = dict()
    for gene, pscores in provean_data.items():
        vobjs = []
        for vobj in uniondata[gene]:
            vchange = vobj.SNPEFF_AMINO_ACID_CHANGE
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
                if in_cobj(vsum.get(v.Gene_ID, []), v) == False:
                    vsum[v.Gene_ID] = vsum.get(v.Gene_ID, []) + [v]
    else:
        for g in variant_list:
            for v in g:
                if in_cobj(vsum.get(v.Gene_ID, []), v) == False:
                    gff = gffdict[v.Gene_ID]
                    v.start = gff.start
                    v.end = gff.end
                    v.strand = gff.strand
                    v.description = gff.description
                    vsum[v.Gene_ID] = vsum.get(v.Gene_ID, []) + [v]
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


def bbh_gene(f, genelist):
    result = dict()
    bdict = bbhdict(f)
    for key in genelist:
        homolog = bdict[key]
        if homolog == "-":
            continue
        else:
            result[homolog] = key
    return result


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
