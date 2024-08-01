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

ls *_dir/*..gz.snpeff.ann.vcf > file.list.txt
head file.list.txt
wc -l file.list.txt
parallel --verbose -j $nt "cp {} $SCRATCHDIR/" :::: file.list.txt # parallelize copying files to the scratch
cd $SCRATCHDIR
echo "files copied to scratch at `date`"
ls
# function that calculates everyting per scaffold and pop parallelized by parallel at the end of the script

count () {
        ploidy=`grep -v "#" $1 | head -n 1 | cut -f 10 | cut -f 1 -d ':' | grep -o "/" | wc -l` # detect ploidy 1 for diploid 3 for tetraploid
        n=$(grep -v "#" $1 | head -n 1 |  tr \\t \\n| wc -l) # get number of individuals
        nind=$(($n-9)) #subtract non-genotype columns

        bcftools filter -e 'ABS(ILEN) > 20' $1 > ${1/..gz./.shortind.} # remove indels longer than 20bp
        bcftools stats -s - ${1/..gz./.shortind.} > ${1/.gz.snpeff.ann.vcf/stats} # extract stats to get sample depth
        grep PSC ${1/.gz.snpeff.ann.vcf/stats} | grep -v "#" | sort -r -k 10 | head -n 6 | cut -f 3 > ${1/.gz.snpeff.ann.vcf/samples} # get list of 6 samples with the highest coverage
        bcftools view -S ${1/.gz.snpeff.ann.vcf/samples} ${1/..gz./.shortind.} | bcftools view --min-ac 1 > ${1/.gz.snpeff.ann.vcf/shortind.subset.ann.vcf} # subset the vcf for the selection of best samples
        vcf=${1/.gz.snpeff.ann.vcf/shortind.subset.ann.vcf}
	bcftools stats -s - $vcf > ${1/.gz.snpeff.ann.vcf/shortind.subset.stats}
        pop=`basename $1 | sed 's/..gz.snpeff.ann.vcf//g'`
        dp_raw=`grep PSC ${1/.gz.snpeff.ann.vcf/stats} | grep -v "#" | sort -r -k 10 | head -n 6 | awk '{s+=$10}END{print s/NR}'` # depth of the input vcf
        dp=`grep PSC ${1/.gz.snpeff.ann.vcf/shortind.subset.stats} | grep -v "#" | sort -r -k 10 | head -n 6 | awk '{s+=$10}END{print s/NR}'` # depth of subset vcf
        sites_tot=`grep -v "#" $vcf | wc -l` # number of sites in subset vcf
        HIGH_tot=`grep -c HIGH $vcf` # 90% of high effect indels are frameshift mutations
        MODERATE_tot=`grep MODERATE $vcf| grep -v HIGH| wc -l` # this category is mainly indels multiple of three
        MODIF_tot=`grep  MODIFIER $vcf| grep -E -v "HIGH|MODERATE|LOW" | wc -l` #mainly indels in non-coding regions
        INTRON_tot=`grep  intron_variant $vcf | grep -E -v "HIGH|MODERATE" | wc -l` #intron indels that are neutral
#       EXON_tot=0 # so far I do not know how to get this number
        frameshift_tot=`grep -c frameshift_variant $vcf`
        scf=`grep -v "#" $vcf | head -n 1 | cut -f 1` # get the scaffold name 


        if [[ $nind -ge 6 ]]
        then
                if [[ $ploidy -eq 1 ]]
                then
                        plo=2 # this part is for diplids
			bgzip $vcf
			tabix -p vcf ${vcf}.gz
                        bcftools view --max-ac 6 ${vcf}.gz > ${vcf/shortind.subset./minor.}  #here consider only sites with minor allele, remove ac>6 it depends on number of indivs you select at line 36
                        vcf=${vcf/shortind.subset./minor.}
                        minor_sites_tot=`grep -v "#" $vcf | wc -l` # add also tot counts for minor dataset
			HIGH_tot_minor=`grep -c HIGH $vcf`
			MODIF_tot_minor=`grep  MODIFIER $vcf| grep -E -v "HIGH|MODERATE|LOW" | wc -l`
			INTRON_tot_minor=`grep  intron_variant $vcf | grep -E -v "HIGH|MODERATE" | wc -l`

                        HIGH_het=`grep HIGH $vcf |cut -f 10-16|grep -E -o "0/1" | wc -l` # get heterozygot genotypes only and count them
                        MODIF_het=`grep MODIFIER $vcf |grep -E -v "HIGH|MODERATE|LOW" |cut -f 10-16|grep -E -o "0/1" | wc -l`
                        HIGH_hom=`grep HIGH $vcf |cut -f 10-16|grep -E -o "1/1" | wc -l` # get homozygous genotypes of given category and count them
                        MODIF_hom=`grep MODIFIER $vcf|grep -E -v "HIGH|MODERATE|LOW"  |cut -f 10-16|grep -E -o "1/1" | wc -l`
                        HIGH_additive=`echo "$HIGH_het+$HIGH_hom*2" |bc` # get number of alleles of given category by summing genotype counts
                        MODIF_additive=`echo "$MODIF_het+$MODIF_hom*2"|bc`
			
			# check also statistics for singletons
			bgzip $vcf
                        tabix -p vcf ${vcf}.gz
			bcftools view --max-ac 1 ${vcf}.gz > ${vcf/minor./singleton.}  #here consider only singletons
			vcf=${vcf/minor./singleton.}
			HIGH_tot_singl=`grep -c HIGH $vcf`
                        MODIF_tot_singl=`grep  MODIFIER $vcf| grep -E -v "HIGH|MODERATE|LOW" | wc -l`


                        printf "%s\t%d\t%s\t%d\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" $pop $plo $scf $nind $dp_raw $dp $sites_tot $HIGH_tot $MODERATE_tot $MODIF_tot $INTRON_tot $frameshift_tot $HIGH_het $HIGH_hom $MODIF_het $MODIF_hom $HIGH_additive $MODIF_additive $minor_sites_tot $HIGH_tot_minor $MODIF_tot_minor $HIGH_tot_singl $MODIF_tot_singl $INTRON_tot_minor

                elif [[ $ploidy -eq 3 ]]
                then
                        plo=4 # this part is for tetraploids
			bgzip $vcf
                        tabix -p vcf ${vcf}.gz
                        bcftools view --max-ac 12 ${vcf}.gz > ${vcf/shortind.subset./minor.} # here consider only sites with minor allele, remove ac>12
                        vcf=${vcf/shortind.subset./minor.}
                        minor_sites_tot=`grep -v "#" $vcf | wc -l`
			HIGH_tot_minor=`grep -c HIGH $vcf`
                        MODIF_tot_minor=`grep  MODIFIER $vcf| grep -E -v "HIGH|MODERATE|LOW" | wc -l`
			INTRON_tot_minor=`grep  intron_variant $vcf | grep -E -v "HIGH|MODERATE" | wc -l`

                        HIGH_het=`grep HIGH $vcf |cut -f 10-16|grep -E -o "0/0/0/1|0/1/1/1|0/0/1/1" | wc -l` #  it counts all dosages
                        HIGH_het1=`grep HIGH $vcf |cut -f 10-16|grep -E -o "0/0/0/1" | wc -l` #  dosage1
                        HIGH_het2=`grep HIGH $vcf |cut -f 10-16|grep -E -o "0/0/1/1" | wc -l` #  dosage2
                        HIGH_het3=`grep HIGH $vcf |cut -f 10-16|grep -E -o "0/1/1/1" | wc -l` #  dosage3
                        MODIF_het=`grep MODIFIER $vcf |grep -E -v "HIGH|MODERATE|LOW" |cut -f 10-16|grep -E -o "0/0/0/1|0/1/1/1|0/0/1/1" | wc -l`
                        MODIF_het1=`grep MODIFIER $vcf |grep -E -v "HIGH|MODERATE|LOW" |cut -f 10-16|grep -E -o "0/0/0/1" | wc -l`
                        MODIF_het2=`grep MODIFIER $vcf |grep -E -v "HIGH|MODERATE|LOW" |cut -f 10-16|grep -E -o "0/0/1/1" | wc -l`
                        MODIF_het3=`grep MODIFIER $vcf |grep -E -v "HIGH|MODERATE|LOW" |cut -f 10-16|grep -E -o "0/1/1/1" | wc -l`
                        HIGH_hom=`grep HIGH $vcf |cut -f 10-16|grep -E -o "1/1/1/1" | wc -l`
                        MODIF_hom=`grep MODIFIER $vcf |grep -E -v "HIGH|MODERATE|LOW" |cut -f 10-16|grep -E -o "1/1/1/1" | wc -l`
                        HIGH_additive=`echo "$HIGH_het1+$HIGH_het2*2+$HIGH_het3*3+$HIGH_hom*4" |bc` # this would probably work for tets
                        MODIF_additive=`echo "$MODIF_het1+$MODIF_het2*2+$MODIF_het3*3+$MODIF_hom*4"|bc`
			bgzip $vcf
                        tabix -p vcf ${vcf}.gz
                        bcftools view --max-ac 1 ${vcf}.gz > ${vcf/minor./singleton.}  #here consider only singletons
                        vcf=${vcf/minor./singleton.}
                        HIGH_tot_singl=`grep -c HIGH $vcf`
                        MODIF_tot_singl=`grep  MODIFIER $vcf| grep -E -v "HIGH|MODERATE|LOW" | wc -l`


                        printf "%s\t%d\t%s\t%d\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" $pop $plo $scf $nind $dp_raw $dp $sites_tot $HIGH_tot $MODERATE_tot $MODIF_tot $INTRON_tot $frameshift_tot $HIGH_het $HIGH_hom $MODIF_het $MODIF_hom $HIGH_additive $MODIF_additive $minor_sites_tot $HIGH_tot_minor $MODIF_tot_minor $HIGH_tot_singl $MODIF_tot_singl $INTRON_tot_minor
                else
                        printf "%s\t%s\n" $pop "strange ploidy"
                fi
        else
                printf "%s\t%s\n" $pop "less than 6 individuals"
        fi
}
export -f count

echo -e "pop\tploidy\tscf\tnind\tdp_raw\tdp\tsites_tot\tHIGH_tot\tMODERATE_tot\tMODIF_tot\tINTRON_tot\tframeshift_tot\tHIGH_het\tHIGH_hom\tMODIF_het\tMODIF_hom\tHIGH_additive\tMODIF_additive\tminor_sites_tot\tHIGH_tot_minor\tMODIF_tot_minor\tHIGH_tot_singl\tMODIF_tot_singl\tINTRON_tot_minor" > genotype.counts.tsv
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
mkdir shortind

mv *minor.ann.vcf.gz minor/
mv *minor.ann.vcf.gz.tbi minor/

mv *shortind.subset.ann.vcf.gz subset/
mv *shortind.subset.ann.vcf.gz.tbi subset/

sfs () { zgrep HIGH $1 |  cut -f 8 | cut -f 1 -d ";" | sed 's/AC=//g' | sort | uniq -c | sort -nr > ${1/vcf.gz/high.sfs}; zgrep MODIFIER $1| grep -E -v "HIGH|MODERATE|LOW" |  cut -f 8 | cut -f 1 -d ";" | sed 's/AC=//g' | sort | uniq -c | sort -nr > ${1/vcf.gz/modifier.sfs}; zgrep intron_variant $1 | grep -E -v "HIGH|MODERATE" | cut -f 8 | cut -f 1 -d ";" | sed 's/AC=//g' | sort | uniq -c | sort -nr > ${1/vcf.gz/intron.sfs} 
}

export -f sfs

echo "constructing sfs `date`"

cd minor # construct sfs from minor vcfs

ls *vcf.gz | parallel --verbose -j $nt "sfs {}"

ls *vcf.gz | cut -f 1 -d "."| sort | uniq > ../pop.txt

# create population level sfs
parallel -j $nt "cat {}*high.sfs > {}.high.sfs" :::: ../pop.txt
parallel -j $nt "cat {}*modifier.sfs > {}.modifier.sfs" :::: ../pop.txt
parallel -j $nt "cat {}*intron.sfs > {}.intron.sfs" :::: ../pop.txt

mkdir sfs
mv *sfs sfs/

cd ../subset # construct sfs from subset vcfs

ls *vcf.gz | parallel -j $nt "sfs {}"

parallel -j $nt "cat {}*high.sfs > {}.high.sfs" :::: ../pop.txt
parallel -j $nt "cat {}*modifier.sfs > {}.modifier.sfs" :::: ../pop.txt
parallel -j $nt "cat {}*intron.sfs > {}.intron.sfs" :::: ../pop.txt


mkdir sfs
mv *sfs sfs/

cd ../

echo "job is done at `date` now move the data back to working directory"

mkdir output.count.hom.het
mkdir singletons.others
mv *vcf* singletons.others/
mv * output.count.hom.het

#tar -cvzf output.count.hom.het.tar.gz output.count.hom.het # ideally you want to parallelize this

find output.count.hom.het -mindepth 1 -maxdepth 1 -type d -print0 | parallel -0 -j $nt tar -czvf {}.tar.gz {} 

mv output.count.hom.het $wd/ || export CLEAN_SCRATCH=false

cd $wd

find output.count.hom.het -type f -name '*.tar.gz' -print0 | parallel -0 -j $nt tar -xzvf {}

echo "finished at `date`"
