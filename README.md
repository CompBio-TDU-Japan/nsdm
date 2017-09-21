# nsdm
NGS Software Development Modules

## Install

```
pip3 install git+https://github.com/CompBio-TDU-Japan/nsdm
```

## function

### nsdm.fileparse

#### vcf_read

return vcf object list.

```
variant_list = nsdm.fileparse.vcf_read(filename)
```

##### vcf object

* alt
* ref
* pos
* annotation
* impact
* gene
* feature

#### gff_read

return gff object dict.

```
variant_list = nsdm.fileparse.gff_read(filename)
```

##### gff object

* start
* end
* strand
* gene
* description
