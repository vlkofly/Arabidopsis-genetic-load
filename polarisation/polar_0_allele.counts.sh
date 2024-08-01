#!/bin/bash
#PBS -N parse_vcf
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=4:mem=60gb:scratch_local=100gb
#PBS -j oe

# This script takes a vcf file and parses it to a format used in program est-sfs for polarisation (acgt format) 
# Several filtering steps precede 

if [ -z "$vcf" ]; then
        echo "Error! Specify relative path (relative to working directory) to vcf "
        echo "qsub -v 'vcf=blablabla.vcf.gz'"
        exit 1
        fi

if [ -z "$out" ]; then
        echo "Specify name of the outfile"
	echo " if not specified default will be used"
	out=${vcf/vcf.gz/allele.count.txt}
        fi

if [ -z "$ref" ]; then
        echo "Specify reference full path, default reference used instead"
        echo "qsub -v '/storage/brno3-cerit/home/filip_kolar/JIC_reference/alygenomes.fasta'"
        ref=/storage/brno3-cerit/home/filip_kolar/JIC_reference/alygenomes.fasta
        fi

if [ -z "$onlymissing" ]; then
        echo "You can run only the missingness command"
        echo "default is to run all the filtering steps"
        onlymissing="no"
        fi


wd=$PBS_O_WORKDIR
nt=4 # number of threads change based on number of scaffolds - but it is better to run one job per a scaffold

refdir=`dirname $ref`
ref=`echo $ref | sed 's/^.*\/\(.*\/.*\)$/\1/g'`

trap 'clean_scratch' TERM EXIT

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 


cp -r $refdir $SCRATCHDIR/
#cp $wd/$vcf $SCRATCHDIR/
cp $vcf $SCRATCHDIR/
#cp $wd/${vcf}.tbi $SCRATCHDIR/ # copy index as well
cp ${vcf}.tbi $SCRATCHDIR/ # copy index as well
#cp $wd/least.missing.individuals.txt $SCRATCHDIR/ # in the case of key generation
vcf=`basename $vcf`


cd $SCRATCHDIR


module add gatk-3.7
module add bcftools-1.6
module add htslib-1.9
tmj="-Djava.io.tmpdir=$SCRATCHDIR"

if [ "$onlymissing" != "yes" ]; then
##select target individuals and only biallelic snps, changed 21.1.2021
java $tmj -XX:ParallelGCThreads=$nt -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R $ref \
       	-V $vcf --excludeNonVariants --restrictAllelesTo BIALLELIC -selectType SNP -o ${vcf/vcf.gz/samplesubs.vcf.gz} -nt $nt -sf least.missing.individuals.txt

##select genotypes with coverage higher than 6
java $tmj -XX:ParallelGCThreads=$nt -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration \
       -R $ref --genotypeFilterExpression 'DP < 6' --genotypeFilterName 'DP' -V ${vcf/vcf.gz/samplesubs.vcf.gz} -o ${vcf/vcf.gz/dp.vcf.gz} -nt $nt

##turn the low coverage genotypes to no-call
java $tmj -XX:ParallelGCThreads=$nt -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration \
      -R $ref --setFilteredGtToNocall -o ${vcf/vcf.gz/dpnc.vcf.gz} -V ${vcf/vcf.gz/dp.vcf.gz} -nt $nt
#remove missingness
java $tmj -XX:ParallelGCThreads=$nt -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants \
	-R $ref --excludeNonVariants --maxNOCALLfraction 0.7 -V ${vcf/vcf.gz/dpnc.vcf.gz} -o ${vcf/vcf.gz/dpnc.m0.7.vcf.gz} -nt $nt # take all sites with 0.5 missingness
else

##in the case you only want to remove missingness

#bgzip -c $vcf > ${vcf}.gz
tabix -p vcf ${vcf}
#vcf=${vcf}.gz
java $tmj -XX:ParallelGCThreads=$nt -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants \
	-R $ref --excludeNonVariants --restrictAllelesTo BIALLELIC -selectType SNP --maxNOCALLfraction 0.7 -V ${vcf} -o ${vcf/vcf.gz/dpnc.m0.7.vcf.gz} -nt $nt # take all sites with 0.5 missingness
fi


### export genotype counts
java $tmj -XX:ParallelGCThreads=$nt -jar $GATK/GenomeAnalysisTK.jar -T VariantsToTable -V ${vcf/vcf.gz/dpnc.m0.7.vcf.gz} \
	-R $ref -F CHROM -F POS -F REF -F ALT -F NSAMPLES -F NCALLED -F AN -F HET -F HOM-REF -F HOM-VAR -GF GT -o ${vcf/vcf.gz/table.tsv}

# turn it into est-sfs format
cp /storage/brno3-cerit/home/filip_kolar/scripts/polar_1_parse.vartbl.py $SCRATCHDIR

# get number of chromosomes
n=`awk '{sum+=$7} END {print sum/(NR-1)}' ${vcf/vcf.gz/table.tsv} |xargs printf '%.0f'`
echo "average number of called alleles: $n "

# subsample to n chromosomes
echo "subsampling the chromosomes to $n"
#python3 polar_1_parse.vartbl.py -tab ${vcf/vcf.gz/table.tsv} -Nallele $n # in the case of populations
python3 polar_1_parse.vartbl.py -tab ${vcf/vcf.gz/table.tsv} -Nallele 56 # in the case to get the polarisation key
echo "number in acgt input:"
cat *acgt.format.tsv | wc -l  

echo "vcf_file numvariants"
for v in *vcf.gz
do
n=`zcat $v | grep -v "#" | wc -l`
echo $v $n
done >> ${vcf/vcf.gz/numbers.ssv}

cat ${vcf/vcf.gz/numbers.ssv}

rm $vcf
#rm *dpnc.vcf*
rm -r JIC_reference  
cp $SCRATCHDIR/* ${wd}/ || export CLEAN_SCRATCH=false

echo "parsing of ${vcf} done"
