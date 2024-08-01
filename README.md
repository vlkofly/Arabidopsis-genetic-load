# Arabidopsis-genetic-load
Scripts used in: Whole genome duplication increases genetic diversity and load in outcrossing Arabidopsis



### adegenet_2023_arenosa.R
Plot PCA and Distance tree using adegenet (Jombart 2008)
### plot.fourfold.div.R
Analyse diversity and Tajima's D (results from [Scantools](https://github.com/mbohutinska/ScanTools_ProtEvol))
### plot.sfs.jan24.R
Plot site frequency spectra (result from Tuomas Hamala [poly_sfs.c script](https://github.com/thamala/polySV) )


## Polarisation using multiple reference genomes
This set of scripts conducts polarisation with program est-sfs (Keightley and Jackson 2018)

### polar_0_allele.counts.sh
Parsing of vcf format
### polar_1_parse.vartbl.py
Generation of acgt format that is input to est-sfs
### polar_2_fix.outgroups.py
Filter alignment of outgroup genomes
### polar_3_est.sfs.sh
Run est-sfs to determine whether allele is ancestral or derived
### plot_polarisation.R
Plot polarisation


## Inference of distribution of fitness effect of new mutations using polyDFE
Population site frequency spectra from est-sfs was used as input to polyDFE (Tataru and Bataillon 2018)

### run_polyDFE.sh
Run several models of polyDFE for each population separately
### polyDFEscript.R
Inspect and plot the results


## Analysis of genetic load
Calculation of load indices for indels, snps and structural variants.

### filter_indels.sh
Filtration of reliable indel variants
### pop.snpeff.sh
Split vcf per population and annotate variants with SNPEff (Cingolani 2022)
### count_hom_het_indels.sh
### count_hom_het_snps.sh
### count_hom_het_indels_norm.sh
### count_hom_het_snps_norm.sh
Count number of sites and homozygous and heterozygous genotypes of given constrain category for indels and snps
### plot.snp.dp.R
### analyse.indels.norm.R
Plot and analyse results of indel load
### analyse.snps.norm.R
Plot and analyse results of snp load
### plot.sv.R
Plot and analyse results of structural variant load



