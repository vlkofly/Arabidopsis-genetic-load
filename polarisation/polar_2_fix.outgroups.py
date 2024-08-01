#!/usr/bin/env python3
# script to fix the outgroup alignment - remove duplicit positions and retain only the one with better score
import argparse
import random
import linecache
import pandas
### Arguments definition:
parser = argparse.ArgumentParser()
parser.add_argument('-bed', type=str, metavar='bed', required=True, help='REQUIRED: bed file to parse; output of MafToBed, without indels')
args = parser.parse_args()

### read through the bedfile 
tab = open(args.bed,'r')
out_name = args.bed.replace(".bed",".fixed.bed")
with open (out_name,'w') as out:
    nt = 0
    total_duplicate=0
    prev=""
    dlines={}
    # maybe collect the duplicit lines into a dictionary
    for l in tab:
        nt += 1
        ls = l.strip("\n")
        #print(ls.split("\t")[0:9])
        CHROM,POS1,POS2,STRAND,SPEC,CHROMS,POSS,BP,STRANDS,SCORE = ls.split("\t")
        line=[CHROM,POS1,POS2,STRAND,SPEC,CHROMS,POSS,BP,STRANDS,SCORE]
        line=list(map(lambda s: s.strip("\t"),line))
        dlines[nt]=line
        #print(line)
        #print(len(dlines))
        if len(dlines) > 1: #compare only when dictionary has two keys
            print(dlines)
            if(dlines[nt-1][1])==(dlines[nt][1]): # if the two consecutive lines are the same position
                print("duplicit lines "+str(dlines[nt-1][1])+" and "+str(dlines[nt][1]))
                #print(str(dlines[nt-1]))
                #print(str(dlines[nt]))
                # take only the lines with the highest qual value but how to understand that there are only two  duplicit lines, or three, and how to take the best out of three?
                # you need to know if the next line is different.
                # the problem is how to get the next line? linecache - come on I could have used it before, but the dictionary solution is neat
                n=linecache.getline(args.bed,nt+1) # it is strangely one based index - get the next line
                # at the end of the file there will be no next position so you need to include if 
                if n == '':
                    print("this is the end of the file")
                    break
                npos=n.split("\t")[1] # extract the chromosomal position of the next line
                #print("next position is "+npos)
                
                if npos != POS1: # if the next position is different (end of ambiguity)
                    
                    v=list(dlines.items())
                    d=pandas.DataFrame.from_dict(dlines,orient="index")
                    d.columns=["CHROM","POS1","POS2","STRAND","SPEC","CHROMS","POSS","BP","STRANDS","SCORE"] # turn the dict into dataframe
                    # select the max value position
                    theline=(d[d.SCORE == d.SCORE.max()]) # get the line with maximal score
                    theline=theline.values.tolist()[0]

                    outline="\t".join(theline)
                    #print (outline)
                    positions=d['POS1']
                    out.write(outline+"\n")
                    print("end of ambiguity with length "+str(len(dlines))+"at positions")
                    print(positions)# write the best aligned position
                    total_duplicate+=1
                    dlines={} # and delete the dictionary



                    #print(d)
                    #


            else: # if the lines are not duplicit delete the dictionary and write the the line to the outfile
                #print("unique lines "+str(dlines[nt-1][1])+" and "+str(dlines[nt][1]))
                outline1="\t".join(dlines[nt-1]) # you need to print both
                outline2="\t".join(dlines[nt]) # you need to print both
                out.write(outline1+"\n")
                # actually you have to remove only the first item in the dictionary to compare all pairs
                #out.write(outline2+"\n")
                del dlines[nt-1] # and remove 
                print(dlines)

print("Fixing of the alignment has finished removing "+str(total_duplicate)+" duplicated sites from "+str(n)+"\n")
