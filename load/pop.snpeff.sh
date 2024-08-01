#!/bin/bash
#PBS -N snpeff_popstats
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=7:mem=100gb:scratch_local=80gb
#PBS -j oe

# This script takes filtered vcf as an argument (vcf) and splits it per population based on popmap which must be supplied as another argument. 
# It annotates variants with SnpEff and stores the resulting vcf into a specified folder and population subfolders.
# supply also name of the snpeff database as a separate argument called db
# Script initially written for mockingbirds but in 2023 modified to use for polyploid plants.

if [ -z "$vcf" ]; then
        echo "Error! Specify vcf file!"
        echo "qsub -v 'vcf=blablabla.vcf.gz'"
        exit 1
        fi

if [ -z "$db" ]; then
        echo "Specify the name of the database in SnpEff!"
        echo "qsub -v 'db=arenosa1'"
        echo "Default will be used: arenosa1"
	db="arenosa1"
        fi

if [ -z "$out" ]; then
        echo "Specify the name of the output directory!"
        echo "qsub -v 'out=arenosa1'"
        echo "default is pop.8x.ac0.indels.annotated"
        out="pop.8x.ac0.indels.annotated"
        fi


if [ -z "$popmap" ]; then
        echo "Error! Specify popmap file!"
	echo "format is POP	sample1,sample2,sampleN"
        echo "qsub -v 'popmap=popgroup.sample.names.all.txt'"
        exit 1
        fi



### define variables ###
nt='6'
wd=$PBS_O_WORKDIR
#popmap="" # file used for splitting vcf format: POP\tsample1,sample2,sampleN, you can use also for splitting individual samples
export outdir=$out # define output folder
# SnpEff directories must contain a database specific for the reference genome used.
export snpeffDB=$db # name of the database
export config="/mnt/storage-brno12-cerit/nfs4/home/filip_kolar/brno3/programs/snpEff.precompiled/snpEff.config"
export bin="/mnt/storage-brno12-cerit/nfs4/home/filip_kolar/brno3/programs/snpEff.precompiled/snpEff.jar"
### add programs ###
module add snpeff/5.1
module add bcftools/1.11
module add htslib/1.9
module add parallel/20160622

trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
cp $wd/$vcf $SCRATCHDIR
cp $wd/${vcf}.tbi $SCRATCHDIR
cp $wd/$popmap $SCRATCHDIR


cd $SCRATCHDIR
echo "files copied to scratch at `date`"
ls

export vcf=`basename $vcf`
popmap=`basename $popmap`

mkdir $outdir

if [[ $vcf != *.gz ]]; then
	bgzip $vcf
	vcf=${vcf/vcf/vcf.gz}
	tabix -p vcf $vcf
fi

export base=${vcf/vcf/}



##function that will be run by parallel per population
popf () {
pop=$1 # first column of popmap
s=$2 # second column of popmap
mkdir ${pop}_dir
echo "processing batch $pop samples $s"
# select samples of a given population 

bcftools view -Oz -s $s --force-samples $vcf | bcftools view -Ov --min-ac 1 -o ${pop}_dir/${pop}.${base}.vcf # consider sites that have at least one non-ref allele in the pop

java -Xmx20g -jar $bin -c $config -s ${pop}_dir/${pop}.${base}.summary.html -csvStats ${pop}_dir/${pop}.${base}.stats.csv -v  $snpeffDB ${pop}_dir/${pop}.${base}.vcf 1> ${pop}_dir/${pop}.${base}.snpeff.ann.vcf 2> ${pop}_dir/${pop}.${base}.snpeff.stderr
}

export -f popf
parallel -j $((nt-1)) --colsep '\t' -a $popmap popf {1} {2} # two columns

# collect the snpeff results from individual populations:
#var="MISSENSE SILENT HIGH MODERATE LOW" # NONSENSE" # collect those variables

#parallel -j 1 "grep {} *_dir/*all.stats.csv |sed 's/\([a-zA-Z]\)\./\1,/g'| sed 's/:/,/g'|sed 's/\//,/g'| tr -d ' ' " ::: $var > pol.snpeff.raw.table.csv # *_dir/*{all,private}.stats.csv replaced because of private missing

#var2="Number_of_variants_processed Number_of_effects Missense_Silent_ratio"
#parallel -j 1 "grep {} *_dir/*all.stats.csv |sed 's/\([a-zA-Z]\)\./\1,/g'| sed 's/:/,/g'|sed 's/\//,/g'| tr -d ' ' " ::: $var2 > pol.snpeff.raw.table2.csv



cp -r * $outdir/
cp -r $outdir $wd/  || export CLEAN_SCRATCH=false
echo " processing of  $vcf done"
