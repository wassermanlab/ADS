##Created by S. Correard
#Updated Feb 9th, 2022

#This script is intended to do the graph of the UMAP colored by the individuals genotypes for a given SNP

library(ggplot2)
library("ggplot2")

current_path_name=getwd()
file_name="1KG_dataset_27022019.GRCh38.phased"
umap.layout_pca_pop=read.table( paste0(current_path_name,file_name,"_umap_layout_pca_pop.txt"), header=TRUE)

###################
scores=read.table(file=paste0("subset_table_score_",file_name,"_allchr"), header = TRUE)

snp_interest_table<- read.table(paste0(current_path_name,"snp_interest"))
snp_interest=snp_interest_table[,1]
scores_snp_interest=which(scores$snp_i==snp_interest)
scores_snp_interest_table=scores[scores_snp_interest,]
write.table(scores_snp_interest_table, file=paste0("variant_interest/scores_snp_interest_table_",snp_interest), quote=FALSE, row.names = FALSE)

Population_score_interest=scores_snp_interest_table[,11]
write.table(Population_score_interest, file="variant_interest/ADS_interest", quote=FALSE, row.names = FALSE)


##Identify the variant of interest

candidates_snp=which(scores$MAF_i>0.257 & scores$MAF_i<0.259 & scores$ADS<0.217)
candidates_snp_table=scores[candidates_snp,]
candidates_snp_table_order=candidates_snp_table[order(candidates_snp_table$ADS,candidates_snp_table$ADS),]
write.table(candidates_snp_table_order, file="candidates_snp", quote=FALSE, row.names = FALSE)


####UMAP representation for a variant of interest

#Load Ped file
ped_file<- read.table(paste0(current_path_name,"snp_interest.ped"), colClasses="character")
colnames(ped_file)=c("FID", "IID", "f", "m", "s", "p", "A1", "A2")
ped_file <- tibble::rowid_to_column(ped_file, "ID")

#Determine the ref and alt allele
#Determine the ref and alt allele
tableA1=table(ped_file$A1)
tableA2=table(ped_file$A2)

#If only one allele (can be 0, 1 or 2) --> next snp
if (nrow(tableA1)==1 & nrow(tableA2)==1){
	next
}

#If 3 columns for both alleles (0, 1 and 2)
if (nrow(tableA1)==nrow(tableA2) & nrow(tableA1)==3){
	table_allele_count=as.data.frame(sort(tableA1+tableA2,decreasing=TRUE))
	Ref_allele=table_allele_count[1,1]
	Alt_allele=table_allele_count[2,1]
	#If Ref allele = 0
	if (Ref_allele==0) {
		Ref_allele=table_allele_count[2,1]
	        Alt_allele=table_allele_count[3,1]
        }
        #If Alt allele = 0
        if (Alt_allele==0) {
                Alt_allele=table_allele_count[3,1]
        }
}

#If 2 columns for both alleles (1 and 2)
if (nrow(tableA1)==nrow(tableA2) & nrow(tableA1)==2){
	table_allele_count=as.data.frame(sort(tableA1+tableA2,decreasing=TRUE))
        Ref_allele=table_allele_count[1,1]
        Alt_allele=table_allele_count[2,1]
        #If only two columns and one of then is 0 (only 0 and 1 or 0 and 2) --> next snp
        if (Ref_allele==0 | Alt_allele==0){
                next
        }
}

#If table 1 has 3 columns (0, 1 and 2) and table 2 has 2 (whichever)
if ((nrow(tableA1)==3) & nrow(tableA1)==(nrow(tableA2)+1)){
	Ref_allele=as.data.frame(sort(tableA1,decreasing=TRUE))[1,1]
        Alt_allele=as.data.frame(sort(tableA1,decreasing=TRUE))[2,1]
        #If Ref allele = 0
        if (Ref_allele==0) {
		Ref_allele=as.data.frame(sort(tableA1,decreasing=TRUE))[2,1]
	        Alt_allele=as.data.frame(sort(tableA1,decreasing=TRUE))[3,1]
        }
        #If Alt allele = 0
        if (Alt_allele==0) {
                Alt_allele=as.data.frame(sort(tableA1,decreasing=TRUE))[3,1]
        }
}

#If table 1 has 2 columns and table 2 has 1
if ((nrow(tableA1)==2) & nrow(tableA1)==(nrow(tableA2)+1)){
	Ref_allele=as.data.frame(sort(tableA1,decreasing=TRUE))[1,1]
	Alt_allele=as.data.frame(sort(tableA1,decreasing=TRUE))[2,1]
        #If only two columns and one of then is 0 (only 0 and 1 or 0 and 2) --> next snp
        if (Ref_allele==0 | Alt_allele==0){
                next
        }
}

#List row for each phenotype
hom1_row=which(ped_file$A1==Ref_allele & ped_file$A2==Ref_allele)
hom2_row=which(ped_file$A1==Alt_allele & ped_file$A2==Alt_allele)
het_row=which(ped_file$A1!=ped_file$A2)
mis_row=which(ped_file$A1==0&ped_file$A2==0)
#het_row=rownames(ped_file[-c(hom1_row, hom2_row, mis_row),])

#List ind for each phenotype
hom1_ind=as.character(ped_file[hom1_row,2])
hom2_ind=as.character(ped_file[hom2_row,2])
het_ind=as.character(ped_file[het_row,2])
mis_ind=as.character(ped_file[mis_row,2])

#Create col next to id with 1:hom ref - 2:het - 3:hom alt - 4:missing
rm("hom1_df", "het_df", "hom2_df", "mis_df")
hom1_df=data.frame(id=hom1_ind, genotype=1)
het_df=data.frame(id=het_ind, genotype=2)

if (length(hom2_ind)>0) {
	hom2_df=data.frame(id=hom2_ind, genotype=3)
}
#mis_df=data.frame(id=mis_ind, genotype=4)
if (exists("mis_df") == TRUE & exists("hom2_df") == TRUE)  {
	ind_geno=rbind(hom1_df,het_df,hom2_df,mis_df)
} else if (exists("mis_df") == FALSE & exists("hom2_df") == TRUE) {
	ind_geno=rbind(hom1_df,het_df,hom2_df)
} else if (exists("mis_df") == TRUE & exists("hom2_df") == FALSE) {
	ind_geno=rbind(hom1_df,het_df,mis_df)
} else if (exists("mis_df") == FALSE & exists("hom2_df") == FALSE) {
	ind_geno=rbind(hom1_df,het_df)
}


#Visualization in UMAP
p12 = ggplot(umap.layout_pca_pop, aes(x_axis_UMAP,y_axis_UMAP)) +
  geom_point(data = umap.layout_pca_pop, aes(x_axis_UMAP,y_axis_UMAP), col = "#FFC107", size = 1, shape=1)+
  geom_point(data = umap.layout_pca_pop[mis_row, ], aes(x_axis_UMAP,y_axis_UMAP), col = "#000000", size = 1, shape=1)+
  geom_point(data = umap.layout_pca_pop[het_row, ], aes(x_axis_UMAP,y_axis_UMAP), col = "#D81B60", size = 1, shape=16)+
  geom_point(data = umap.layout_pca_pop[hom2_row, ], aes(x_axis_UMAP,y_axis_UMAP), col = "#1E88E5", size = 1, shape=16)+
  theme_void() +
  labs(title = paste0( "UMAP - ",snp_interest))

#Print the representation
ggsave(filename = paste0("variant_interest/UMAP_",snp_interest,".png"), p12, dpi =300)

##Graph UMAP score distribution for given MAF
MAF_score_interest=scores_snp_interest_table[,10]

#min MAF value
min_MAF=MAF_score_interest-(MAF_score_interest/10)

#Max MAF value
max_MAF=MAF_score_interest+(MAF_score_interest/10)

#Subset table for values in this MAF
scores_MAF_interest=which(scores$MAF_i>min_MAF & scores$MAF_i<max_MAF)
scores_MAF_interest_table=scores[scores_MAF_interest,]
write.table(scores_MAF_interest_table, file="scores_MAF_interest_table", quote=FALSE, row.names = FALSE)

#Create the graph
p13=ggplot(scores_MAF_interest_table, aes(ADS)) +
  geom_histogram(binwidth=.0002, colour="black", fill="white") +
  theme_minimal()+
#  abline(v=c(Population_score_interest), col=c("red"))+
  geom_vline(xintercept = Population_score_interest , linetype = "dashed", colour = "red")+
  labs(title =paste0( "ADS distribution - ",min_MAF," - ",max_MAF, ",", snp_interest))

#Print the representation
ggsave(filename = paste0("variant_interest/score_distrib_MAF_",snp_interest,".png"), p13, dpi =300)










