#!/usr/bin/env python3
import nsdm
def main():
    a = nsdm.fileparse.vcf_read("/Users/yuto/Dropbox/Lab/code/ngsdata/allvcf/Generation01.vcf")
    print(a)

if __name__ == '__main__':
    main()
