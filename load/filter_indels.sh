#!/bin/bash
#PBS -N filter_indels
#PBS -l walltime=40:00:00
#PBS -l select=1:ncpus=2:mem=150gb:scratch_local=100gb
#PBS -j oe

# This script takes raw variants after joint genotyping as a an argument (vcf) together with a relative path to a reference genome (ref).
# It exports filtered indels. There are two filtering steps, first filtration is based on GATK best practices with excess depth as well
# Second filtering is based on number of reads supporting each indel within a genotype.
# Specific samples are also selected, defined in sample.txt file.

if [ -z "$vcf" ]; then
	echo "Error! Specify relative path (relative to working directory) to vcf file!"
        echo "qsub -v 'vcf=blablabla.vcf.gz'"
        exit 1
        fi


if [ -z "$ref" ]; then
        echo "Error! Specify relative path (relative to working directory) to reference fasta!"
        echo "qsub -v 'ref=arenosa1.fna'"
        exit 1
        fi



wd=$PBS_O_WORKDIR

trap 'clean_scratch' TERM EXIT

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 


cp $wd/$vcf $SCRATCHDIR/
cp $wd/${vcf}.tbi $SCRATCHDIR
cp $wd/${ref} $SCRATCHDIR
cp $wd/${ref}.fai $SCRATCHDIR
cp $wd/samples.txt $SCRATCHDIR # sample file, one smaple a line
ref=`basename $ref`
vcf=`basename $vcf`


cd $SCRATCHDIR


module add bcftools/1.11
module add htslib/1.9

# subset only the individuals that are used in scantools
# filtering gatk best practices plus excess depth 2x mean depth
# get the average depth
#dp=`bcftools view -v indels -S samples.txt --force-samples $vcf | bcftools query -f '%INFO/DP\n' | awk '{s+=$1} END {print s/NR}'` # it is enough to calculat depth only for a subset of sites
#dp=`echo "$dp*2" | bc` # times two
#echo "excess depth level is $dp"

# -a does not work due to some PL incompatibility
# multiallelic sites are split into separate lines
#bcftools view -v indels --force-samples -S samples.txt $vcf | bcftools norm -f $ref --force -m-|bcftools view -v indels | tee ${vcf/vcf.gz/raw.norm.indels.vcf}| bcftools filter -e "FS>40.0 || SOR>3 || MQ<40 || MQRankSum<-5.0 || MQRankSum>5.0 || QD<2.0 || ReadPosRankSum<-4.0 || INFO/DP > $dp" > ${vcf/vcf.gz/filtered.indels.vcf}

#bcftools filter -S . -i 'FORMAT/AD[*:1] >= 8' ${vcf/vcf.gz/filtered.indels.vcf} | bcftools filter -e 'AC=0' > ${vcf/vcf.gz/filtered.8x.indels.vcf}

#use this for filter 5x dataset
bcftools filter -S . -i 'FORMAT/AD[*:1] >= 5' ${vcf} | bcftools filter -e 'AC=0' > ${vcf/vcf/filtered.5x.indels.vcf}

bcftools stats -d8,58,40 -F $ref ${vcf/vcf.gz/filtered.5x.indels.vcf} > ${vcf/vcf.gz/filtered.5x.indels.stats}


echo "vcf_file numvariants"
for v in `ls *vcf`
do
n=`cat $v | grep -v "#" | wc -l`
echo $v $n
done >> ${vcf/vcf/5x.numbers.ssv}

#bgzip ${vcf/filtered.indels/filtered.5x.indels}


rm $vcf

cp $SCRATCHDIR/* ${wd}/ || export CLEAN_SCRATCH=false

#cp $SCRATCHDIR/* $wd/filtervcf_toxo/ || export CLEAN_SCRATCH=false
echo "filtering of ${vcf} done"
