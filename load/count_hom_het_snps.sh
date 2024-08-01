#!/bin/bash
#PBS -N count_hom_het
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=45:mem=100gb:scratch_local=400gb
#PBS -j oe

# this script goes to population directories that are created by pop.snpeff.sh and it calculates numbers of sites/genotypes of specific fitness category
# it works with diploids and tetraploids for vcf split per scaffold


nt='44'
wd=${PBS_O_WORKDIR}

module rm parallel/20160622
module add bcftools/1.11
module add htslib/1.9
module add parallel/20200322 

echo "started at `date`"

cd $wd  #go to working directory 
echo "Working directory is:"
pwd
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi

cp /storage/brno12-cerit/home/filip_kolar/brno3/scripts/load/plot.snp.dp.R $SCRATCHDIR/

ls *_dir/*..gz.snpeff.ann.vcf > file.list.txt
head file.list.txt
wc -l file.list.txt
parallel --verbose -j $nt "cp {} $SCRATCHDIR/" :::: file.list.txt # parallelize copying files to the scratch
cd $SCRATCHDIR
echo "files copied to scratch at `date`"
ls
# function that calculates everyting per scaffold and pop parallelized by parallel at the end of the script

count () {
        pop=`echo $1 | cut -f 1 -d "." `
	scf=`grep -v "#" $1 | head -n 1 | cut -f 1` # get the scaffold name 
        ploidy=`grep -v "#" $1 | head -n 1 | cut -f 10 | cut -f 1 -d ':' | grep -o "/" | wc -l` # detect ploidy 1 for diploid 3 for tetraploid
        n=$(grep -v "#" $1 | head -n 1 |  tr \\t \\n| wc -l) # get number of individuals
        nind=$(($n-9)) #subtract non-genotype columns
	#bcftools query -f '[%DP \n]' $1 | sort | uniq -c > DP.distr.${pop}.${scf}.txt # check the distribution of genotype depth
	#cat DP.distr.${pop}.${scf}.txt |  sed 's/^[[:space:]]*//g'| sponge DP.distr.${pop}.${scf}.txt
	#Rscript plot.snp.dp.R DP.distr.${pop}.${scf}.txt $pop $scf >> rscript.log


        bcftools filter -S . -e 'FMT/DP<8' $1 | bcftools view -e 'F_MISSING > 0.5' --min-ac 1 > ${1/..gz./.8x.} # what additional filtering should I use for snps, yes missingness 0.5
	# I was thinking about QD but maybe I will use again some hard filter 8x 
        bcftools stats -s - ${1/..gz./.8x.} > ${1/.gz.snpeff.ann.vcf/stats} # extract stats to get sample depth
        grep PSC ${1/.gz.snpeff.ann.vcf/stats} | grep -v "#" | sort -r -k 10 | head -n 6 | cut -f 3 > ${1/.gz.snpeff.ann.vcf/samples} # get list of 6 samples with the highest coverage
        bcftools view -S ${1/.gz.snpeff.ann.vcf/samples} ${1/..gz./.8x.} | bcftools view --min-ac 1 > ${1/.gz.snpeff.ann.vcf/8x.subset.ann.vcf} # subset the vcf for the selection of best samples
        vcf=${1/.gz.snpeff.ann.vcf/8x.subset.ann.vcf}

	bcftools stats -s - $vcf > ${1/.gz.snpeff.ann.vcf/8x.subset.stats}
        dp_raw=`grep PSC ${1/.gz.snpeff.ann.vcf/stats} | grep -v "#" | sort -r -k 10 | head -n 6 | awk '{s+=$10}END{print s/NR}'` # depth of the input vcf
        dp=`grep PSC ${1/.gz.snpeff.ann.vcf/8x.subset.stats} | grep -v "#" | sort -r -k 10 | head -n 6 | awk '{s+=$10}END{print s/NR}'` # depth of subset vcf

	#count different categories of snps
        sites_tot=`grep -v "#" $vcf | wc -l` # number of sites in subset vcf
        MISSENSE_tot=`grep -c missense_variant $vcf` # Non-synonymous mutations, based on SNPEff their effect is moderate 
	NONSENSE_tot=`grep stop_gained $vcf| wc -l` # stop_gained
        SILENT_tot=`grep  synonymous_variant $vcf| grep -E -v "HIGH|MODERATE" | wc -l` # Synonymous mutations, based on SNPeff the effec is LOW 
        WEAK_tot=`grep  MODIFIER $vcf | grep -E -v "HIGH|MODERATE|LOW" | wc -l` #Weak category is mainly intron variants and non-cds variants
        scf=`grep -v "#" $vcf | head -n 1 | cut -f 1` # get the scaffold name 


        if [[ $nind -ge 6 ]]
        then
                if [[ $ploidy -eq 1 ]]
                then
                        plo=2 # this part is for diplids
			bgzip $vcf
			tabix -p vcf ${vcf}.gz
                        bcftools view --max-ac 6 ${vcf}.gz > ${vcf/8x.subset./minor.}  #here consider only sites with minor allele, remove ac>6 it depends on number of indivs you select at line 36
                        vcf=${vcf/8x.subset./minor.}
                        minor_sites_tot=`grep -v "#" $vcf | wc -l` #counts for minor dataset
			MISSENSE_tot_minor=`grep -c missense_variant $vcf`
			NONSENSE_tot_minor=`grep stop_gained $vcf| wc -l` 
			SILENT_tot_minor=`grep  synonymous_variant $vcf| grep -E -v "HIGH|MODERATE" | wc -l`
			WEAK_tot_minor=`grep  MODIFIER $vcf | grep -E -v "HIGH|MODERATE|LOW" | wc -l`

                        MISSENSE_het=`grep missense_variant $vcf |cut -f 10-16|grep -E -o "0/1" | wc -l` # get heterozygot genotypes only and count them
                        SILENT_het=`grep synonymous_variant $vcf |grep -E -v "HIGH|MODERATE" |cut -f 10-16|grep -E -o "0/1" | wc -l`
                        MISSENSE_hom=`grep missense_variant $vcf |cut -f 10-16|grep -E -o "1/1" | wc -l` # get homozygous genotypes of given category and count them
                        SILENT_hom=`grep synonymous_variant $vcf|grep -E -v "HIGH|MODERATE"  |cut -f 10-16|grep -E -o "1/1" | wc -l`
                        MISSENSE_additive=`echo "$MISSENSE_het+$MISSENSE_hom*2" |bc` # get number of alleles of given category by summing genotype counts
                        SILENT_additive=`echo "$SILENT_het+$SILENT_hom*2"|bc`

			#NONSENSE_het= I will try just missense silent 
			#NONSENSE_hom=
			#NONSENSE_additive=


			
			# check also statistics for singletons
			bgzip $vcf
                        tabix -p vcf ${vcf}.gz
			bcftools view --max-ac 1 ${vcf}.gz > ${vcf/minor./singleton.}  #here consider only singletons
			vcf=${vcf/minor./singleton.}
			MISSENSE_tot_singl=`grep -c missense_variant $vcf`
                        SILENT_tot_singl=`grep  synonymous_variant $vcf| grep -E -v "HIGH|MODERATE" | wc -l`


                        printf "%s\t%d\t%s\t%d\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" $pop $plo $scf $nind $dp_raw $dp $sites_tot $MISSENSE_tot $NONSENSE_tot $SILENT_tot $WEAK_tot $MISSENSE_het $MISSENSE_hom $SILENT_het $SILENT_hom $MISSENSE_additive $SILENT_additive $minor_sites_tot $MISSENSE_tot_minor $SILENT_tot_minor $NONSENSE_tot_minor $MISSENSE_tot_singl $SILENT_tot_singl $WEAK_tot_minor

                elif [[ $ploidy -eq 3 ]]
                then
                        plo=4 # this part is for tetraploids
			bgzip $vcf
                        tabix -p vcf ${vcf}.gz
                        bcftools view --max-ac 12 ${vcf}.gz > ${vcf/8x.subset./minor.} # here consider only sites with minor allele, remove ac>12
                        vcf=${vcf/8x.subset./minor.}
                        minor_sites_tot=`grep -v "#" $vcf | wc -l`
			MISSENSE_tot_minor=`grep -c missense_variant $vcf`
			NONSENSE_tot_minor=`grep stop_gained $vcf| wc -l` 
                        SILENT_tot_minor=`grep  synonymous_variant $vcf| grep -E -v "HIGH|MODERATE" | wc -l`
			WEAK_tot_minor=`grep  MODIFIER $vcf | grep -E -v "MISSENSE|NONSENSE" | wc -l`

                        MISSENSE_het=`grep missense_variant $vcf |cut -f 10-16|grep -E -o "0/0/0/1|0/1/1/1|0/0/1/1" | wc -l` #  it counts all dosages
                        MISSENSE_het1=`grep missense_variant $vcf |cut -f 10-16|grep -E -o "0/0/0/1" | wc -l` #  dosage1
                        MISSENSE_het2=`grep missense_variant $vcf |cut -f 10-16|grep -E -o "0/0/1/1" | wc -l` #  dosage2
                        MISSENSE_het3=`grep missense_variant $vcf |cut -f 10-16|grep -E -o "0/1/1/1" | wc -l` #  dosage3
                        SILENT_het=`grep synonymous_variant $vcf |grep -E -v "HIGH|MODERATE" |cut -f 10-16|grep -E -o "0/0/0/1|0/1/1/1|0/0/1/1" | wc -l`
                        SILENT_het1=`grep synonymous_variant $vcf |grep -E -v "HIGH|MODERATE" |cut -f 10-16|grep -E -o "0/0/0/1" | wc -l`
                        SILENT_het2=`grep synonymous_variant $vcf |grep -E -v "HIGH|MODERATE" |cut -f 10-16|grep -E -o "0/0/1/1" | wc -l`
                        SILENT_het3=`grep synonymous_variant $vcf |grep -E -v "HIGH|MODERATE" |cut -f 10-16|grep -E -o "0/1/1/1" | wc -l`
                        MISSENSE_hom=`grep missense_variant $vcf |cut -f 10-16|grep -E -o "1/1/1/1" | wc -l`
                        SILENT_hom=`grep synonymous_variant $vcf |grep -E -v "HIGH|MODERATE" |cut -f 10-16|grep -E -o "1/1/1/1" | wc -l`
                        MISSENSE_additive=`echo "$MISSENSE_het1+$MISSENSE_het2*2+$MISSENSE_het3*3+$MISSENSE_hom*4" |bc` # this would probably work for tets
                        SILENT_additive=`echo "$SILENT_het1+$SILENT_het2*2+$SILENT_het3*3+$SILENT_hom*4"|bc`
			bgzip $vcf
                        tabix -p vcf ${vcf}.gz
                        bcftools view --max-ac 1 ${vcf}.gz > ${vcf/minor./singleton.}  #here consider only singletons
                        vcf=${vcf/minor./singleton.}
                        MISSENSE_tot_singl=`grep -c missense_variant $vcf`
                        SILENT_tot_singl=`grep  synonymous_variant $vcf| grep -E -v "HIGH|MODERATE" | wc -l`


                        printf "%s\t%d\t%s\t%d\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" $pop $plo $scf $nind $dp_raw $dp $sites_tot $MISSENSE_tot $NONSENSE_tot $SILENT_tot $WEAK_tot $MISSENSE_het $MISSENSE_hom $SILENT_het $SILENT_hom $MISSENSE_additive $SILENT_additive $minor_sites_tot $MISSENSE_tot_minor $SILENT_tot_minor $NONSENSE_tot_minor $MISSENSE_tot_singl $SILENT_tot_singl $WEAK_tot_minor
                else
                        printf "%s\t%s\n" $pop "strange ploidy"
                fi
        else
                printf "%s\t%s\n" $pop "less than 6 individuals"
        fi
}
export -f count

echo -e "pop\tploidy\tscf\tnind\tdp_raw\tdp\tsites_tot\tMISSENSE_tot\tNONSENSE_tot\tSILENT_tot\tWEAK_tot\tMISSENSE_het\tMISSENSE_hom\tSILENT_het\tSILENT_hom\tMISSENSE_additive\tSILENT_additive\tminor_sites_tot\tMISSENSE_tot_minor\tSILENT_tot_minor\tNONSENSE_tot_minor\tMISSENSE_tot_singl\tSILENT_tot_singl\tWEAK_tot_minor" > genotype.counts.tsv
#ls */*.gz.snpeff.ann.vcf|head -n 1 | par.gz.l --bar -j 9 "count {}" >> genotype.counts.tsv
#ls */*snpeff.ann.vcf | par.gz.l --bar -j 9 "count {}" >> genotype.counts.tsv
#ls BAB_dir/*..gz.snpeff.ann.vcf | parallel --bar -j 14 "count {}" >> genotype.counts.BAB.tsv
ls *vcf | parallel  -j $nt "count {}" >> genotype.counts.tsv

mv genotype.counts.tsv $wd/ || export CLEAN_SCRATCH=false # genotype counts goes back to the working directory
# copy also the filtered vcfs back to the working directory for further inspection
rm *..gz.snpeff.ann.vcf

# construct sfs from minor vcf
echo "constructing sfs started at `date`"
mkdir minor
mkdir subset 

mv *minor.ann.vcf.gz minor/
mv *minor.ann.vcf.gz.tbi minor/

mv *8x.subset.ann.vcf.gz subset/
mv *8x.subset.ann.vcf.gz.tbi subset/

sfs () { zgrep missense_variant $1 |  cut -f 8 | cut -f 1 -d ";" | sed 's/AC=//g' | sort | uniq -c | sort -nr > ${1/vcf.gz/missense.sfs}; zgrep synonymous_variant $1| grep -E -v "HIGH|MODERATE" |  cut -f 8 | cut -f 1 -d ";" | sed 's/AC=//g' | sort | uniq -c | sort -nr > ${1/vcf.gz/silent.sfs} 
}

export -f sfs

echo "constructing sfs `date`"

cd minor # construct sfs from minor vcfs

ls *vcf.gz | parallel --verbose -j $nt "sfs {}"

ls *vcf.gz | cut -f 1 -d "."| sort | uniq > ../pop.txt

# create population level sfs and save them in specific folder
mkdir pop.sfs
parallel -j $nt "cat {}*missense.sfs > pop.sfs/{}.missense.sfs" :::: ../pop.txt
parallel -j $nt "cat {}*silent.sfs > pop.sfs/{}.silent.sfs" :::: ../pop.txt

mkdir sfs
mv *sfs sfs/

cd ../subset # construct sfs from subset vcfs

ls *vcf.gz | parallel -j $nt "sfs {}"
mkdir pop.sfs
parallel -j $nt "cat {}*missense.sfs > pop.sfs/{}.missense.sfs" :::: ../pop.txt
parallel -j $nt "cat {}*silent.sfs > pop.sfs/{}.silent.sfs" :::: ../pop.txt

mkdir sfs
mv *sfs sfs/

cd ../

echo "job is done at `date` now move the data back to working directory"

mkdir output.count.hom.het
mkdir singletons.others
mv *scaffold* singletons.others/
mv * output.count.hom.het

#tar -cvzf output.count.hom.het.tar.gz output.count.hom.het # ideally you want to parallelize this

mkdir $wd/output.count.hom.het

find output.count.hom.het -mindepth 1 -maxdepth 1 -type d -print0 | parallel -0 -j $nt "tar -cvf {}.tar {} && mv {}.tar $wd/output.count.hom.het" 

mv -u output.count.hom.het $wd/
export CLEAN_SCRATCH=false

cd $wd

find output.count.hom.het -type f -name '*.tar' -print0 | parallel -0 -j $nt "tar -xvf {} && rm {}"

echo "finished at `date`"
