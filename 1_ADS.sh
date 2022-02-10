#!/bin/sh
#SBATCH -o outfile_ADS
#SBATCH -e errfile_ADS
#SBATCH --time=27-0
#SBATCH --mem=40G
#SBATCH --ntasks-per-node=4

#Created by Solenne Correard
#Updated on Feb 9th, 2022

#This pipeline does several things, in order :
#1. Prepare the data for the UMAP using plink: Select randomly a susbset of variant on each chromosome and merge them into one file that will be used to calculate the PCA and create the UMAP (particulary useful if the vcf are splitted by chromosomes or if there is too many variants to calcule le PCA quickly).
#2. Prepare the data for the calculation of the ADS using plink: Keep only the data with AC>1, and split in small subset of 250,000 variants to allow parallelization and faster processing by R.
#3. Generate the UMAP using R (Refer to the R script for more details)
#4. Calcule the ADS for each variant using R  (Refer to the R script for more details)


module load StdEnv/2018.3 plink/1.9b_5.2-x86_64 nixpkgs/16.09  gcc/7.3.0 r/4.0.2

#Files and path 
#Path to the input file
input_dir="/PATH/TO/INPUT/FILES"
#Output name : name wanted for the output data
output_file="1KG_dataset_27022019.GRCh38.phased"

#Name of the input file (if the dataset is splited by chr, replace the chr value by ${i} in the name
#file_name="your_file_name.vcf.gz"
#To change in the loop if the input dataset is splited by chr

#Parameters to define
#thin_threshold is the probability of retaining each sample at random. A thin_threshold of 0.05 (default value) keeps randomly 5% of the variants.
thin_threshold="0.05"
#number_PCs is the number of principal components that plink must calculate based on the thined data. The PCs will be used to generate the UMAP (default value : 15)
number_PCs="15"

mkdir -p splited_by_chr_thin
mkdir -p splited_by_chr_AC2
mkdir -p splited_by_chr_AC2/full_chr
mkdir -p graph_output

for i in {1..22}
do
	#File name if the input dataset is splited by chr
	file_name="ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
	#Thin the data
	plink --thin ${thin_threshold} --set-missing-var-ids @:#[b37]\$1,\$2 --vcf ${input_dir}/${file_name} --recode --out splited_by_chr_thin/${output_file}_thin${thin_threshold}_chr$i
	
	#Keep the varaints with AC>1 and split the chr in smaller pieces for UMAP calculation
	#Plink : Keep all the variant with AC>1 by chr
	plink --min-ac 2 --chr $i --set-missing-var-ids @:#[b37]\$1,\$2 --vcf ${input_dir}/${file_name} --recode --out splited_by_chr_AC2/${output_file}_AC2_chr$i
	#Count the number of variants with AC>1 by chr
	wc -l splited_by_chr_AC2/${output_file}_AC2_chr$i.map > splited_by_chr_AC2/wc-l_chr$i
	nb_snps=$(awk '{print $1;}' splited_by_chr_AC2/wc-l_chr$i)

	echo "chr"
	echo $i
	echo ${nb_snps}

	k=1
	while [[ $k -le ${nb_snps} ]]
	do
		echo "k"
		echo $k
	       	l=$((k+249999))
		m=$((l+1))
		echo "l"
		echo $l

		#Write the list of variants to keep in each sub dataset
		sed -n "${k},${l}p;${m}q" splited_by_chr_AC2/${output_file}_AC2_chr${i}.map > splited_by_chr_AC2/snps_chr${i}_${k}_${l}
        	
		#Extract the variants from the list previously created
		plink --chr ${i} --extract splited_by_chr_AC2/snps_chr${i}_${k}_${l} --file splited_by_chr_AC2/${output_file}_AC2_chr${i} --recode --out splited_by_chr_AC2/${output_file}_AC2_chr${i}_${k}_${l}
	
		((k = k + 250000))
	done
done

#Merge the thined data to create the PCA and the UMAP
#Create a file with all the ped/map names from chr2
for j in {2..22}
do
	echo splited_by_chr_thin/${output_file}_thin${thin_threshold}_chr$j.ped	splited_by_chr_thin/${output_file}_thin${thin_threshold}_chr$j.map >> splited_by_chr_thin/list_ped_map_thin${thin_threshold}
done

#Merge the file with thined data from chr 1 with the list of file created previously
plink --file splited_by_chr_thin/${output_file}_thin${thin_threshold}_chr1 --merge-list splited_by_chr_thin/list_ped_map_thin${thin_threshold} --recode 12 --out ${output_file}_thin${thin_threshold}

#Calculate the Principal components that will be used to generate the UMAP
plink --file ${output_file}_thin${thin_threshold} --pca ${number_PCs} --out ${output_file}_thin${thin_threshold}_PCA

#Move the full schoromosomes ped and map files in another folder
mv splited_by_chr_AC2/1KG_dataset_27022019.GRCh38.phased_AC2_chr?.* splited_by_chr_AC2/full_chr/
mv splited_by_chr_AC2/1KG_dataset_27022019.GRCh38.phased_AC2_chr??.* splited_by_chr_AC2/full_chr/

#Start the UMAP generation and the ADS calculation using a R script
Rscript 1_R_ADS.R




