#!/usr/bin/env python3

import argparse
import random
### Arguments definition:
parser = argparse.ArgumentParser()
parser.add_argument('-tab', type=str, metavar='tab', required=True, help='REQUIRED: table from VariantsToTable GATK command ')
parser.add_argument('-Nallele',type=int,metavar='Nallele',required=False,default=100, help='Specify number of alleles to which the complete dataset will be subsampled. Default is 100')

args = parser.parse_args()
### generate list of sites that needs to be polarised

### read through the table 
tab = open(args.tab,'r')

out_name = args.tab.replace(".table.tsv",".acgt.format.tsv")
with open (out_name,'w') as out:
    nt = 0
    for l in tab:
        if l.startswith("CHROM"):
            print("skipping headline")
        else:
            nt += 1
            ls = l.strip("\n")
            #print(ls.split("\t")[0:9])
            CHROM,POS,REF,ALT,NSAMPLES,NCALLED,HET,HOMREF,HOMVAR = ls.split("\t")[0:9]
            GT=list(ls.split("\t")[10:])
            gtall="_".join(GT)
            print(GT)
            ###you need to randomly sample 100 alleles from the dataset and avoid missing genotypes
            allist=[] # generate list of alleles per site
            for g in GT:
                gg=g.split("/")
                for a in gg:
                    if a != ".": # remove missing alleles
                        allist.append(a)
            nchr=len(allist)
            print("number of chroms:"+str(nchr)+" number to subsample to "+str(args.Nallele))
            # make the random selection of alleles
            if nchr >= args.Nallele:
                print("enough chromosomes to select from:"+str(nchr))
                randomlist=random.sample(allist,args.Nallele)
                randomstring="_".join(randomlist)
            else:
                print("not enough chroms to select from:"+str(nchr))
                continue
            nchrom=0
            if len(REF) > 1:
                print("skipping insertion at "+POS)
                continue
            out.write(CHROM+"\t"+str(int(POS)-1)+"\t"+POS+"\t")

            for b in "A" "C" "G" "T":
                ac=0 # now count number of each allele per site
                ac=randomstring.count(b)
                nchrom+=ac # how many chromosomes sampled per site?
                if b != "T":
                    out.write(str(ac)+',')
                else:
                    out.write(str(ac)+'\n')
                    print(POS,nchrom) 


#
