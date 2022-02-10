#!/bin/sh
#SBATCH -o outfile_variant_interest
#SBATCH -e errfile_variant_interest
#SBATCH --time=03:45:00
#SBATCH --mem=40G
#SBATCH --ntasks-per-node=4

#Created by Solenne Correard
#Updated on Feb 9th, 2022

#This pipeline allow to draw the UMAP while coloring each individual depending on their genotype for a given SNP

module load StdEnv/2018.3 plink/1.9b_5.2-x86_64

output_file="1KG_dataset_27022019.GRCh38.phased"

mkdir -p variant_interest

##Identify in which ped/map the snp of interest is located
#Low MAF, Low score
#MAF = 0.001177394 , score = 0.002737324

#snp_interest="22:36265860[b37]A,G"

#Or load several SNPs from a file
#22:36265860[b37]A,G
#9:93709790[b37]C,G

while IFS= read -r line
do
	snp_interest=$line	

	pos_interest=$(echo "$snp_interest" | cut -f 1 -d '[')
	echo $snp_interest > snp_interest
	file_snp_interest=$(grep -l $pos_interest splited_by_chr_AC2/${output_file}_AC2_chr*.map)
	
	echo $pos_interest
	echo $snp_interest
	echo $file_snp_interest
	#Remove ".map" of the file_name
	file_snp_interest_short=$(echo $file_snp_interest | rev | cut -f 2- -d '.' | rev)
	echo $file_snp_interest_short

	#Extract the data for the snp of interest
	module load StdEnv/2018.3 plink/1.9b_5.2-x86_64
	plink --file $file_snp_interest_short -extract snp_interest --recode --out snp_interest

	#Create image
	module load StdEnv/2020 r/4.0.2

	Rscript R_3_variant_interest.R

done < list_snp_interest

#Merge all the table score in one table score
awk 'FNR==1 && NR!=1{next;}{print}' variant_interest/scores_snp_interest_table_* > variant_interest/scores_snps_interest

#Cleaning
rm snp_interest*





