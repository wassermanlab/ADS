#!/bin/sh
#SBATCH -o outfile_annotation
#SBATCH -e errfile_annotation
#SBATCH --time=01:45:00
#SBATCH --mem=40G
#SBATCH --ntasks-per-node=4

#Created by Solenne Correard
#Updated on Feb 9th, 2022

#This pipeline annotate the variants for which a ADS was calculated to plot the score in function of the annotation
module load StdEnv/2018.3 nixpkgs/16.09  gcc/7.3.0 r/4.0.2 annovar

#annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene humandb/

#Files and path 
#Name of the input file (file obtained from script 1_ADS.sh with the different scores
file_name="1KG_dataset_27022019.GRCh38.phased"
input_file=subset_table_score_${file_name}_allchr

#Split second column to have a annovar like format
##Need chr pos pos Ref Alt
awk -F '[][,:" "]' '{print $2 "\t" $3 "\t" $3 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18}' ${input_file} > ${input_file}_annovar_input

annotate_variation.pl --outfile ${input_file}_annotated --geneanno -buildver hg38 ${input_file}_annovar_input humandb/

rm ${input_file}_annovar_input

#Start the graph generation using a R script
Rscript R_2_annotate_variants.R
