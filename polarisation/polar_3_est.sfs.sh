#!/bin/bash
#PBS -N polarise 
#PBS -l select=1:ncpus=1:mem=100gb:scratch_local=100gb
#PBS -l walltime=10:00:00
#PBS -j oe

# supply argument chr=Chr1 & config & outd 
# polarise per chromosome, because to put whole genome on one pile would be too heavy
# example of running the script qsub -v 'chr=Chr1' ~/scripts/polarisation3_est_sfs.sh


# the script has been modified for the purpose of arabidopsis polarisation 29.7.2020


# you have to modify the line 80 -100 how to parse the outgroups
module add gsl-2.1-gcc
module add bedtools-2.26.0
module add python34-modules-gcc
export CPATH=/software/gsl/2.1/include/:$CPATH
export LIBRARY_PATH=/software/gsl/2.1/lib/:$LIBRARY_PATH

wd=$PBS_O_WORKDIR
bin="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/programs/est-sfs-release-2.03/est-sfs"
beddir="/storage/brno3-cerit/home/filip_kolar/600_genomes/polarisation/outgr_alignments" # the output from tba pipeline
# for file definition use relative paths 

if [ -z "$outd" ]; then
        echo "outidr not specified usig default est_sfs_results"
	outd="est_sfs_results"
fi

if [ -z "$config" ]; then
        echo "configuration file not specified usig default config_est_sfs_kimura.txt"
	config="config_kimura.txt"
fi

if [ -z "$focfreq" ]; then
       echo "You need to supply focal freq file if you do not want default:"
       focfreq=$wd/focal_freq_sfs.txt
fi

if [ -z "$outgr" ]; then
       echo "You need to supply bedfile with outgroups"
       outgr=$wd/focal_freq_sfs.txt
       exit 1
fi


if [ -z "$chr" ]; then
	echo "You need to supply chrname"
	exit 1
	fi

# keep the chr as well 



# the vcf was processed this way to match est-sfs format which is:
# 20,0,0,0        0,0,0,1 0,0,0,1 0,0,0,1 # ACGT counts

trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi

cp $wd/$focfreq $SCRATCHDIR/
cp $wd/$outgr $SCRATCHDIR/
cp $wd/$config $SCRATCHDIR/ # est_sfs config file should be in the working directory
cp $wd/seedfile.txt $SCRATCHDIR/ # the same as seedfile.txt
cd $SCRATCHDIR

#bed=`ls *bed`
focfreq=`basename $focfreq`
outgr=`basename $outgr`


# remove indels in the bed alignment when sorting do not forget to redirect tmp
awk '{if ($8 ~ /.,.,./) print $0}' $outgr | sort -T $SCRATCHDIR -k1,1 -k2,2n > ${outgr/.bed/.noindel.bed}

# fix the outgroup alignment (ambiguities)
cp ~/scripts/polar_2_fix.outgroups.py .
python3 polar_2_fix.outgroups.py -bed ${outgr/.bed/.noindel.bed} > fix.${outgr/.bed/}.log

outgr=${outgr/.bed/.noindel.fixed.bed} # the script renames the file


# let's assume the both beds are sorted and do the intersection when you want to keep all sites that are in vcf
# because est-sfs allows missingness in outgroups

intersectBed -sorted -loj -a $focfreq -b $outgr > ${chr}.intersected

# parse the intersected file so it conforms to the format of est-sfs
# get the column with outgroup values corresponds to melanotis,ficedula,sturnus

### this option for outgroup melanotis,sturnus,ficedula
#paste <(cut -f 4 ${chr}.intersected) <(cut -f 12 ${chr}.intersected | awk -F, '{am=0;tm=0;gm=0;cm=0;as=0;ts=0;gs=0;cs=0;af=0;tf=0;gf=0;cf=0; \
#if($1 ~ /[aA]/)  am=1;if($1 ~ /[tT]/) tm=1;if($1 ~ /[gG]/) gm=1;if($1 ~ /[cC]/) cm=1; if($3 ~ /[aA]/) as=1; if($3 ~ /[tT]/) ts=1; \
#if($3 ~ /[gG]/) gs=1;if($3 ~ /[cC]/) cs=1;if($2 ~ /[aA]/) af=1;if($2 ~ /[tT]/) tf=1;if($2 ~ /[gG]/) gf=1;if($2 ~ /[cC]/) cf=1; \
#printf "%d,%d,%d,%d %d,%d,%d,%d %d,%d,%d,%d\n", am,cm,gm,tm,as,cs,gs,ts,af,cf,gf,tf }') > ${chr}.est_sfs.input
#missing=`grep "0,0,0,0 0,0,0,0 0,0,0,0" ${chr}.est_sfs.input | wc -l`

### this option for thaliana and capsella, capsella is first:
paste <(cut -f 4 ${chr}.intersected) <(cut -f 12 ${chr}.intersected | awk -F, '{am=0;tm=0;gm=0;cm=0;as=0;ts=0;gs=0;cs=0;af=0;tf=0;gf=0;cf=0; \
if($1 ~ /[aA]/)  am=1;if($1 ~ /[tT]/) tm=1;if($1 ~ /[gG]/) gm=1;if($1 ~ /[cC]/) cm=1; if($3 ~ /[aA]/) as=1; if($3 ~ /[tT]/) ts=1; \
if($3 ~ /[gG]/) gs=1;if($3 ~ /[cC]/) cs=1;if($2 ~ /[aA]/) af=1;if($2 ~ /[tT]/) tf=1;if($2 ~ /[gG]/) gf=1;if($2 ~ /[cC]/) cf=1; \
printf "%d,%d,%d,%d %d,%d,%d,%d\n", as,cs,gs,ts,af,cf,gf,tf }') > ${chr}.est_sfs.input

### this option for polyglottos,sturnus,ficedula
# one more column was added
#paste <(cut -f 4,5 ${chr}.intersected) <(cut -f 13 ${chr}.intersected | awk -F, '{am=0;tm=0;gm=0;cm=0;as=0;ts=0;gs=0;cs=0;af=0;tf=0;gf=0;cf=0; \
#if($1 ~ /[aA]/)  am=1;if($1 ~ /[tT]/) tm=1;if($1 ~ /[gG]/) gm=1;if($1 ~ /[cC]/) cm=1; if($3 ~ /[aA]/) as=1; if($3 ~ /[tT]/) ts=1; \
#if($3 ~ /[gG]/) gs=1;if($3 ~ /[cC]/) cs=1;if($2 ~ /[aA]/) af=1;if($2 ~ /[tT]/) tf=1;if($2 ~ /[gG]/) gf=1;if($2 ~ /[cC]/) cf=1; \
#printf "%d,%d,%d,%d %d,%d,%d,%d\n", as,cs,gs,ts,af,cf,gf,tf }') > ${chr}.est_sfs.input


# how many sites are not covered by any outgroup?
missing=`grep "0,0,0,0 0,0,0,0" ${chr}.est_sfs.input | wc -l`

# now run the est-sfs
$bin $config ${chr}.est_sfs.input seedfile.txt ${chr}.sfs.output.txt ${chr}.site.output.txt


# get sites to be switched the same as polarisation key that I have for lyrata
# check the distribution of probabilities and based on that derive thresholds
# there will be some sites that will be clear adepts for repolarisation, approx when the prob (column 3 in sfs output) of major allele being ancestral will be <0.4
# than some gray zone in betwee 0.4-0.6 sites that cannot be clearly polarised and then >0.6 sites where major is ancestral (that will be majority of sites)

#paste <(cut -f 1-3,13 ${chr}.intersected) ${chr}.est_sfs.input <(tail -n +9 ${chr}.site.output.txt) |  tr -s " " \\t  > ${chr}.site.complete.output.bed #to assign the prob estimates with chrom positions
paste <(cut -f 1-3,12 ${chr}.intersected) ${chr}.est_sfs.input <(tail -n +9 ${chr}.site.output.txt) |  tr -s " " \\t  > ${chr}.site.complete.output.bed # in the case without polyglottos as outgroup

# add header in the Rscript


# check also if you change the outgroups what will happen with sfs and polarisation key

cat ${chr}.sfs.output.txt | tr "," \\n | sponge ${chr}.sfs.output.txt

# now write an R script that would summarise the polarisation and pick the polarisation key
Rscript ~/scripts/plot_polarisation.R ${chr}.site.complete.output.bed ${chr}.sfs.output.txt

# I will need to annotate alleles in vcf which is derived and which ancestral.


wc -l *
echo "while there is $missing sites without reference"
#rm *bed
#rm focal_freq*
mkdir $outd
mv * $outd
cp -r $outd $wd/ || export CLEAN_SCRATCH=false

