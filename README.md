# Allele Dispersion Score (ADS)

## Goal of the ADS

The Allele Dispersion Score (ADS) was created to address the challenge of  describing the dispersion of an allele across populations

For each allele, the ADS can range from 0 (the variant is limited to a subset of close individuals within the whole cohort) to 1 (the variant is spread among the individuals present in the cohort).

The ADS is associated with each allele, complementing the overall allele frequency. As for the allele frequency, the ADS will differ depending on the cohort studied.



## Calculation of the ADS

The ADS is calculated based on projected distance between individuals carrying the variant within a Uniform Manifold Approximation and Projection (UMAP).

The calculation of the ADS is a two steps process: (i) a UMAP is created based on a random subset of genotype data from all individuals. (ii) for each variant, the ADS is calculated based on the projected distance between the individuals carrying the alternate allele. 

You can refer to the manuscript for more details.



## Data and code availability

### To calculate the ADS with your own dataset

To generate the ADS for your own dataset, you can use the two script available in this folder :
1_ADS.sh (including PLINK steps) and R_1_ADS.R (including UMAP generation and the ADS calculation)

The ADS pipeline can be used with vcf, vcf.gz or ped/map files as input files.

The current scripts are optimized to use on a server with a SLURM manager. A Nextflow version of the ADS pipeline should be available in the coming months.



### To reproduce the mauscript results

The dataset used to create and test the ADS is available on the IGSR website : [1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/).

The code used to calculate the ADS on this dataset and compare to the label score is available in the "IGSR_dataset_Manuscript/code/" folder. It is composed of a shell script (1_ADS.sh, including PLINK steps) and an R script (R_1_ADS.R, including UMAP generation and the ADS calculation). It includes the figure generation and the correlation calculation.

The code used to create the other figures of the manuscript are available in two additional couple of scripts :

- 2_annotate_variants.sh and R_2_annotate_variants.R : Annotate the variants and plot the ADS per annotation (Figure 2D)

- 3_snp_interest.sh and R_3_variant_interest.R : For a list of SNP of interest, it generates the UMAP with individuals colored according to their genotype for the SNP of interest and the ADS distribution for variants with a comparable MAF.

Each script is highly commented to explain the role of each step, but feel free to contact us if you have any question or issue with reproducing the results.
If you regenerate the ADS with the same dataset from scratch, you are likely going to obtain slightly different results as the UMAP will be slightly different. However, the overall results and conslusion should be the same. 

These scripts were processed on  [Compute Canada](https://www.computecanada.ca) using a SLURM scheduler.



### Manuscript results

The files resulting from the ADS calculation on the IGSR dataset are available in the "IGSR_dataset_Manuscript/results/" folder. Some intermediate files were too big to be included in the repo, please, contact us if you need them. 

The final dataset with the ADS for each variant is splitted by chr, under the name "subset_table_score_1KG_dataset_27022019.GRCh38.phased_chr{chr_numer}".




## Contact

For any question or comment about the tool, feel free to contact the first author of the manuscript : Solenne Correard : scorreard@cmmt.ubc.ca

## Citation

If you are generating the ADS for you dataset, please, cite our manuscript :
[To add]
