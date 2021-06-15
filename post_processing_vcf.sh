
#Reheader the vcf and change the samples names :
bcftools query -l FB_HC_CommonVars_Merged_116F2_DP20_AO10_genome_13B.vcf.gz > old_samples.txt

#edit with awk the "old_samples.txt" to new sample names: 

bcftools reheader -s rename.samples.txt -o Renamed_SC-F2.vcf.gz --threads 4 FB_HC_CommonVars_Merged_116F2_DP20_AO10_genome_13B.vcf.gz

#bcftools to norm the vcf file and remove any duplicated markers: 

bcftools norm -d all --threads=4 Renamed_SC-F2.vcf.gz -O z -o Renamed_Norm_SC-F2.vcf.gz

#Beagle to phase snps statisitcally:

beagle gt=Renamed_Norm_SC-F2.vcf.gz out=Renamed_Phased_SC-F2.vcf.gz burnin=5 iterations=15 ne=250000 nthreads=4

#Vcftools final filtering out of indels and phenotyped clones:

vcftools --gzvcf Renamed_Phased_SC-F2.vcf.gz --remove-indels --keep phenotypes.txt --recode --stdout |bgzip > Renames_Norm_Phase_Filtered.vcf.gz

#plink to LD prune to ~10,000 SNPs:

plink --vcf Renames_Norm_Phase_Filtered.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 1000 500 0.88 --out Renames_Norm_Phase_Filtered
#PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
#(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
#Logging to Renames_Norm_Phase_Filtered.log.
#Options in effect:
#  --allow-extra-chr
#  --double-id
#  --indep-pairwise 1000 500 0.88
#  --out Renames_Norm_Phase_Filtered
#  --set-missing-var-ids @:#
#  --vcf Renames_Norm_Phase_Filtered.vcf.gz

#15990 MB RAM detected; reserving 7995 MB for main workspace.
#--vcf: Renames_Norm_Phase_Filtered-temporary.bed +
#Renames_Norm_Phase_Filtered-temporary.bim +
#Renames_Norm_Phase_Filtered-temporary.fam written.
#353853 variants loaded from .bim file.
#353853 missing IDs set.
#110 people (0 males, 0 females, 110 ambiguous) loaded from .fam.
#Ambiguous sex IDs written to Renames_Norm_Phase_Filtered.nosex .
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 110 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.999934.
#353853 variants and 110 people pass filters and QC.
#Note: No phenotypes present.
#Pruned 86511 variants from chromosome 1, leaving 3125.
#Pruned 38551 variants from chromosome 2, leaving 1327.
#Pruned 14154 variants from chromosome 3, leaving 418.
#Pruned 32488 variants from chromosome 4, leaving 952.
#Pruned 17203 variants from chromosome 5, leaving 287.
#Pruned 39407 variants from chromosome 6, leaving 1162.
#Pruned 38022 variants from chromosome 7, leaving 1003.
#Pruned 14336 variants from chromosome 8, leaving 214.
#Pruned 10614 variants from chromosome 9, leaving 359.
#Pruned 24294 variants from chromosome 10, leaving 373.
#Pruned 15357 variants from chromosome 11, leaving 325.
#Pruned 12837 variants from chromosome 12, leaving 534.
#Pruning complete.  343774 of 353853 variants removed.
#Marker lists written to Renames_Norm_Phase_Filtered.prune.in and
#Renames_Norm_Phase_Filtered.prune.out .

#total: 10,079 SNPs


#Turn 0/1 vcf genotypes into AB genotypes for R/qtl

zcat Renames_Norm_Phase_Filtered_Pruned.vcf.gz |grep -v "#"|cut -f 10-120 |sed -e "s/\///g"|sed -e 's/00/AA/g'|sed -e 's/11/BB/g'|sed -e 's/01/AB/g'|sed -e 's/10/BA/g' > QTL_genotypes

#get maker names:

awk -F ":" '{print"M_"$1"_"$2}' Renames_Norm_Phase_Filtered.prune.in > marker_names

#put together:

paste marker_names QTL_genotypes > Rqtl_input.temp

#Add the individuals into the first column in Excel:

#Make the other files needed:

#phyiscal map:

cat marker_names |while read line; do echo $line |tr "\n" "\t"; echo $line|awk -F "R" '{print$2}'|awk -F "_" '{print$1"\t"$2}'; done > SC_physical.temp

#gentic map: 

cat marker_names |while read line; do echo $line |awk -F "_" '{print$2"\t"$3}' > TEMP; cat TEMP|while read item; do grep -w "$item" marker.locations |sed -e "s/CHR//g";done; done > SC_genetic.temp


paste marker_names SC_genetic.temp |awk -F "\t" '{print$1"\t"$2"\t"$4"\t"$5}' > temp && mv temp SC_genetic.temp

#edit both in excel to be the correct format

