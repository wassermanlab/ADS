##Created by S. Correard
#Updated Feb 9th, 2022

#This script is intended to calculate the Allele Dispersion Score (ADS) based on a UMAP for all the variants in a file
#Previous script : Plink --> Load a vcf or ped/map files and create ped/map files splited by chr
		#Plink --> Create pca

#Current script : Generate the UMAP based on the PCs previously calculated and calculate the ADS based on UMAP
#Steps
#2 R --> Create UMAP [Combine the n dimensions in 2 dim] for distance score calculation and vizualisation (need pca data of step 1)
#3 R --> Create distance matrix (distance between each ind) for distance score calculation and vizualisation (need UMAP of step 2)
#4 R --> Calcule the ADS for each SNP (need the matrix distance from 3, and ped/map splited by chr) 
#5 R --> Combine the tables for each chr to have one final table


#files
current_path_name=getwd()
file_name="1KG_dataset_27022019.GRCh38.phased"

#umap fonction
library(umap)
#ggplot function
library(ggplot2)
#melt function
library(reshape)
#top_n function
library(dplyr)
library(svglite)
library(cowplot)

#2 R --> Create UMAP [Combine the n dimensions in 2 dim] for distance score calculation and vizualisation (need PCA data from previous step)

#OR load UMAP if created previously
#umap.layout=read.table("umap_layout.txt", header=TRUE)
#umap.layout_pca=read.table( paste0(current_path_name,file_name,"_umap_layout_pca.txt"), header=TRUE)
#umap.layout_pca_pop=read.table( paste0(current_path_name,file_name,"_umap_layout_pca_pop.txt"), header=TRUE)

pca_data<- read.table(paste0(current_path_name,file_name,"_thin0.05_PCA.eigenvec"), header = FALSE)


###UMAP
#Umap configuration
custom.config = umap.defaults
custom.config$metric="euclidean" #character or function; determines how distances between data points are computed. When using a string, available metrics are: euclidean, manhattan. Other available generalized metrics are: cosine, pearson, pearson2. Note the triangle inequality may not be satisfied by some generalized metrics, hence knn search may not be optimal. When using metric.function as a function, the signature must be function(matrix, origin, target) and should compute a distance between the origin column and the target columns  (Correard et al., 2021 : euclidian)
custom.config$min_dist=0.5 #numeric; determines how close points appear in the final layout (Correard et al., 2021 : 0.5)
custom.config$n_neighbors=15 #integer; number of nearest neighbors (Correard et al., 2021 : 15)
custom.config$n_components=2 #integer; dimension of target (output) space (Correard et al., 2021 : 2)
custom.config$input="data" #character, use either "data" or "dist"; determines whether the primary input argument to umap() is treated as a data matrix or as a distance matrix (Correard et al., 2021 : data)

#Select the first 15 PCs to do the UMAP
data.umap_pca=pca_data[,c(3:17)]

#Calculate the UMAP
umap_pca=umap(data.umap_pca, custom.config)
umap.layout_pca=as.data.frame(umap_pca$layout)
colnames(umap.layout_pca)=c("x_axis_UMAP", "y_axis_UMAP")
write.table(umap.layout_pca, file=paste0(current_path_name,file_name,"_umap_layout_pca.txt"), quote=FALSE, row.names = FALSE)
umap.labels = pca_data[, "V1"]
umap.layout_labels_pca=cbind(umap.layout_pca, umap.labels)


#Vizualisation of the UMAP
p=ggplot(umap.layout_labels_pca, aes(x_axis_UMAP,y_axis_UMAP)) +
  theme_void() +
  geom_point(size = 0.3) +
  labs(title =paste0( "UMAP, ", file_name))

#Print the representation
ggsave(filename = paste0(current_path_name,file_name,"_UMAP.svg"), p, dpi =300 )



#3 R --> Create distance matrix for UMAP (distance between each ind) for distance score calculation and vizualisation (need UMAP of step 2)
#Matrice des distance
dmat_UMAP  <-  as.matrix(stats::dist(umap.layout_pca, diag = FALSE, upper = FALSE))
dist_UMAP.melt=dist(umap.layout_pca, diag = FALSE, upper = FALSE)

#Order the matrix values in increasing order to obtain max and min sums
dist_UMAP.melt_sorted_min_first=sort(dist_UMAP.melt)
dist_UMAP.melt_sorted_max_first=sort(dist_UMAP.melt, decreasing = TRUE)

write.table(dist_UMAP.melt_sorted_max_first, paste0(current_path_name,file_name,"_dist_UMAP.melt_sorted_max_first.txt"), quote=FALSE, row.names = FALSE)


#4 R --> Calcule the ADS for each SNP (need the nmatrix distance from 3, and ped/map splited by chr) 
#Info
#ind = Individual 
#Hom1 = homozygous Reference
#Hom2 = Homozygous alternatif
#Het = Heterozygous
#Mis = missing (0/0)


#Loop for each chr
map_file=list.files(path=paste0(current_path_name, "splited_by_chr_AC2/"), pattern=paste0(file_name,"_AC2_chr.+.map"))

#Loop for each chr
for (l in 1:length(map_file)){
	
	#Clean objects
	rm(Alt_allele_i,col_ped_A1,col_ped_A2,dist_i,dist_ind_het,dist_ind_het_hom2,dist_ind_hom2,dmat.ind_het,dmat.ind_het_hom2,dmat.ind_hom2,het_row_i,hom1_row_i,hom2_row_i,ind_het_i,ind_hom1_i,ind_hom2_i,ind_mis_i,k_het_hom2_i,k_het_i,k_hom2_i,MAF_i,max_dist_i,min_dist_i,mis_row_i,n_allele_alt,n_allele_ref,n_allele_tot,n_ind_het_hom2,n_ind_het_i,n_ind_hom1_i,n_ind_hom2_i,n_ind_mis_i,norm_UMAP_score,pos_i,Ref_allele_i,snp_i,subset_ped_snpi,table_allele_count,tableA1,tableA2) 


	#Load Map/ped files
	#map_file_l<- read.table(paste0(current_path_name,"splited_by_chr_AC2/",file_name_AC,"_chr",l,".map"))
	map_file_l<-read.table(paste0(current_path_name,"splited_by_chr_AC2/",map_file[l]))
	colnames(map_file_l)=c("chr", "snp", "pos(mg)", "pos(bp)")
	file_name_for_ped=substr(map_file[l],1,nchar(map_file[l])-4)
	ped_file_l<- read.table(paste0(current_path_name,"splited_by_chr_AC2/",file_name_for_ped,".ped"), colClasses="character")
  
	##Loop to split the number of variant
	slots=c(seq(0, nrow(map_file_l)-1, by=50000), nrow(map_file_l))
  
	#For each slot of variants
	for (j in 1:(length(slots)-1)){
		min_i=slots[j]+1
		max_i=slots[j+1]
  
		## Loop to calculate the distance score for each variant
		table_dist_UMAP_score=data.frame()
	
		#for (i in 1:nrow(map_file_l)){
		for (i in min_i:max_i){
			snp_i=as.character(map_file_l[i,2])
			chr_i=map_file_l[i,1]
			pos_i=map_file_l[i,4]
			col_ped_A1=6+2*i-1
			col_ped_A2=col_ped_A1+1
			subset_ped_snpi=ped_file_l[,c(1,2,col_ped_A1,col_ped_A2)]
			colnames(subset_ped_snpi)=c("FID", "IID", "A1", "A2")
    
			#Determine the ref and alt allele
			tableA1=table(subset_ped_snpi$A1)
			tableA2=table(subset_ped_snpi$A2)
	
			#If only one allele (can be 0, 1 or 2) --> next snp
			if (nrow(tableA1)==1 & nrow(tableA2)==1){
				next
			} 
      
			#If 3 columns for both alleles (0, 1 and 2)
			if (nrow(tableA1)==nrow(tableA2) & nrow(tableA1)==3){
				table_allele_count=as.data.frame(sort(tableA1+tableA2,decreasing=TRUE))
				Ref_allele_i=table_allele_count[1,1]
				Alt_allele_i=table_allele_count[2,1]
				#If Ref allele = 0
				if (Ref_allele_i==0) {
					Ref_allele_i=table_allele_count[2,1]
					Alt_allele_i=table_allele_count[3,1]
				}
				#If Alt allele = 0
				if (Alt_allele_i==0) {
					Alt_allele_i=table_allele_count[3,1]
				}
			}	
      
			#If 2 columns for both alleles (1 and 2)
			if (nrow(tableA1)==nrow(tableA2) & nrow(tableA1)==2){
        			table_allele_count=as.data.frame(sort(tableA1+tableA2,decreasing=TRUE))
        			Ref_allele_i=table_allele_count[1,1]
        			Alt_allele_i=table_allele_count[2,1]
        			#If only two columns and one of then is 0 (only 0 and 1 or 0 and 2) --> next snp
        			if (Ref_allele_i==0 | Alt_allele_i ==0){
          				next
        			}	
			}
      
			#If table 1 has 3 columns (0, 1 and 2) and table 2 has 2 (whichever)
			if ((nrow(tableA1)==3) & nrow(tableA1)==(nrow(tableA2)+1)){
        			Ref_allele_i=as.data.frame(sort(tableA1,decreasing=TRUE))[1,1]
        			Alt_allele_i=as.data.frame(sort(tableA1,decreasing=TRUE))[2,1]
        			#If Ref allele = 0
        			if (Ref_allele_i==0) {
        				Ref_allele_i=as.data.frame(sort(tableA1,decreasing=TRUE))[2,1]
					Alt_allele_i=as.data.frame(sort(tableA1,decreasing=TRUE))[3,1]
        			}	
        			#If Alt allele = 0
        			if (Alt_allele_i==0) {
        				Alt_allele_i=as.data.frame(sort(tableA1,decreasing=TRUE))[3,1]
        			}
			}

			#If table 1 has 2 columns and table 2 has 1
 			if ((nrow(tableA1)==2) & nrow(tableA1)==(nrow(tableA2)+1)){
        			Ref_allele_i=as.data.frame(sort(tableA1,decreasing=TRUE))[1,1]
        			Alt_allele_i=as.data.frame(sort(tableA1,decreasing=TRUE))[2,1]
        			#If only two columns and one of then is 0 (only 0 and 1 or 0 and 2) --> next snp
        			if (Ref_allele_i==0 | Alt_allele_i ==0){
					next
				}
			}
      
			#List ind / row for each genotype
			hom1_row_i=which(subset_ped_snpi$A1==Ref_allele_i & subset_ped_snpi$A2==Ref_allele_i)
			hom2_row_i=which(subset_ped_snpi$A1==Alt_allele_i & subset_ped_snpi$A2==Alt_allele_i)
			het_row_i=which(subset_ped_snpi$A1!=subset_ped_snpi$A2)
			mis_row_i=which(subset_ped_snpi$A1==0&subset_ped_snpi$A2==0)
	      
			#Define each category
			ind_hom1_i=c(hom1_row_i)
			ind_hom2_i=c(hom2_row_i)
			ind_het_i=c(het_row_i)
			ind_mis_i=c(mis_row_i)
      	
			#Number of ind in each category
			n_ind_hom1_i=length(ind_hom1_i)
			n_ind_hom2_i=length(ind_hom2_i)
			n_ind_het_i=length(ind_het_i)
			n_ind_mis_i=length(mis_row_i)
			n_ind_het_hom2=n_ind_het_i+n_ind_hom2_i
			n_allele_ref=2*length(ind_hom1_i)+length(ind_het_i)
			n_allele_alt=2*length(ind_hom2_i)+length(ind_het_i)
			n_allele_tot=n_allele_ref+n_allele_alt
			MAF_i=n_allele_alt/n_allele_tot
      	
			#If no individual or only one individual carry the variant, then no score calculated, variant is skipped
			if (n_ind_het_hom2==0 | n_ind_het_hom2==1){
	        		next
	        	} else if (MAF_i>0.5){
	        		next
	        	} else {
	
				#Calculate the ADS for the allele of interest
	        		#Depend on the genotype of the individual
        	
				#Subset distance matrix (with only individuals hom alt for the variant of interest)
        			dmat.ind_hom2=dmat_UMAP[(ind_hom2_i),c(ind_hom2_i)]
        			dist_ind_hom2=2*sum(dmat.ind_hom2)
				#number of cells considered
				k_hom2_i=(n_ind_hom2_i^2-n_ind_hom2_i)/2 
      	
        			#Subset distance matrix (with only individuals het for the variant of interest)
        			dmat.ind_het=dmat_UMAP[(ind_het_i),c(ind_het_i)]
        			dist_ind_het=sum(dmat.ind_het)/2
        			#number of cells considered
        			k_het_i=(((n_ind_het_i^2)-n_ind_het_i)/2)
      	
        			#Subset distance matrix (with only individuals hom alt in one side of the matrix and individuals het on the other side)
        			dmat.ind_het_hom2=dmat_UMAP[(ind_hom2_i),c(ind_het_i)]
        			dist_ind_het_hom2=2*sum(dmat.ind_het_hom2)
        			#number of cells considered
        			k_het_hom2_i=n_ind_hom2_i*n_ind_het_i
     	 
        			#Add the values with allele coefficient
        			dist_i=dist_ind_hom2+dist_ind_het_hom2+dist_ind_het
      
        			##Min value (For normalization)
        			min_dist_i=4*sum(dist_UMAP.melt_sorted_min_first[(1:k_hom2_i)])+2*sum(dist_UMAP.melt_sorted_min_first[((k_hom2_i+1):(k_hom2_i+1+k_het_hom2_i))])+sum(dist_UMAP.melt_sorted_min_first[((k_hom2_i+1+k_het_hom2_i+1):(k_hom2_i+1+k_het_hom2_i+1+k_het_i))])

	 			##Max value (For normalization)
				max_dist_i=4*sum(dist_UMAP.melt_sorted_max_first[(1:k_hom2_i)])+2*sum(dist_UMAP.melt_sorted_max_first[((k_hom2_i+1):(k_hom2_i+1+k_het_hom2_i))])+sum(dist_UMAP.melt_sorted_max_first[((k_hom2_i+1+k_het_hom2_i+1):(k_hom2_i+1+k_het_hom2_i+1+k_het_i))])

        
        			#Normalized data
        			ADS = (dist_i-min_dist_i)/(max_dist_i-min_dist_i)
			}
        
			#Write table with snp of interest coordinates, genotype distribution, MAF and UMAP score
			#chr, snp, pos, ref allele, alt allele, num of ind/cat (hom1, hom2, het, mis), UMAP score
			temp.table_dist_UMAP_score=cbind(chr_i, snp_i, pos_i, as.character(Ref_allele_i), as.character(Alt_allele_i), n_ind_hom1_i, n_ind_hom2_i, n_ind_het_i, n_ind_mis_i, MAF_i, ADS,  dist_i, max_dist_i, min_dist_i)
			table_dist_UMAP_score=rbind.data.frame(table_dist_UMAP_score,temp.table_dist_UMAP_score)
		}
  
		#Write table slots
		write.table(table_dist_UMAP_score, file=paste0("table_score_",file_name_for_ped,"_slot", j), quote=FALSE, row.names = FALSE)
	}
  
	#Write table 
	list_tables_slots <- list.files(pattern = paste0("table_score_",file_name_for_ped, "_slot"))
	tables_dist_score_slots=lapply(list_tables_slots, read.table, header=TRUE)
	combined_tables_dist_score_slots=do.call(rbind,tables_dist_score_slots)
	write.table(combined_tables_dist_score_slots, file=paste0("table_score_",file_name_for_ped), quote=FALSE, row.names = FALSE)
	file.remove(list_tables_slots)
}


#6 R --> Combine the tables for each chr to have one final table
list_tables_dist_score <- list.files(pattern = paste0("table_score_",file_name,"_AC2_chr"))
tables_dist_score=lapply(list_tables_dist_score, read.table, header=TRUE)
combined_tables_dist_score=do.call(rbind,tables_dist_score)
combined_tables_dist_score_ordered= combined_tables_dist_score[order(combined_tables_dist_score$chr_i, combined_tables_dist_score$pos_i),]
write.table(combined_tables_dist_score_ordered, file=paste0("table_score_",file_name,"_allchr"), quote=FALSE, row.names = FALSE)
#file.remove(list_tables_dist_score)

##Subset the final tab
combined_tables_dist_score_ordered=read.table(file=paste0("table_score_",file_name,"_allchr"), header=TRUE)
col_to_keep=c("chr_i", "snp_i", "pos_i", "V4", "V5", "n_ind_hom1_i", "n_ind_hom2_i", "n_ind_het_i", "n_ind_mis_i", "MAF_i", "ADS")
subset_table= combined_tables_dist_score_ordered[col_to_keep]
write.table(subset_table, file=paste0("subset_table_score_",file_name,"_allchr"), quote=FALSE, row.names = FALSE)


##Figures
#scores=read.table(file=paste0("table_score_",file_name,"_allchr"), header = TRUE)
scores=subset_table
options(scipen = 999)

#calculated MAF distribution
p1=ggplot(scores, aes(MAF_i)) +
  geom_histogram(binwidth=.0002, colour="black", fill="grey") +
  theme_minimal()+
  labs(title = "Calculated MAF distribution")

ggsave(filename = paste0("graph_output/",file_name,"_MAF_distribution.svg"), p1, dpi =320 )

#ADS distribution
p6=ggplot(scores, aes(ADS)) +
  theme_minimal()+
  geom_histogram(binwidth=.01, colour="grey75", fill="grey90") +
  geom_density(aes(y=0.01 * ..count..), colour="black") +
  labs(title = "ADS distribution", x= "Allele Dispersion Score (ADS)", y=" ")+
  theme(text = element_text(family = 'Helvetica Neue'),
        plot.title = element_text(size = 19),
        axis.title = element_text(size = 19),
        axis.text = element_text(size = 15))

ggsave(filename = paste0("graph_output/",file_name,"_ADS_distribution.svg"), p6, dpi =320 )


#Boxplot for distribution of the ADS per MAF
#Create the bins
scores_bin = scores %>% mutate( bin=cut_width(MAF_i, width=0.01, boundary=0) )
#Calcule the number of observation per bin
data_per_bin <- count(scores_bin, bin)

#Do the plot
p22=ggplot(scores_bin, aes(x=bin, y=ADS)) +
        geom_boxplot()+
        theme_minimal()+
	theme(text = element_text(family = 'Helvetica Neue'),
        	plot.title = element_text(size = 19),
        	axis.title = element_text(size = 19),
        	axis.text = element_text(size = 15),
        	panel.grid.major.x = element_blank() ,
        	panel.grid.major.y = element_line( size=.1, color="grey85"))+
#	theme(axis.text.x  = element_text(angle=90, vjust=0.5))+
	labs( title ="ADS distribution", subtitle = "MAF 0 to 0.5", x="Minor Allele Frequency (MAF)", y="Allele Dispersion Score (ADS)")

#Too big to save as svg
ggsave(filename = paste0("graph_output/",file_name,"_ADS_distrib_allMAF.png"), p22, dpi= 300 )

print("Number of variants per bin for supp data")
print(data_per_bin)

#Boxplot for distribution of the ADS for low MAF
#Select low MAF
scores_MAF_interest=which( scores$MAF_i<0.05)
scores_MAF_interest_table=scores[scores_MAF_interest,]
#Create the bins
scores_small_MAF_bin = scores_MAF_interest_table %>% mutate( bin=cut_width(MAF_i, width=0.001) )

#Calcule the number of observation per bin
data_per_bin_small_MAF <- count(scores_small_MAF_bin, bin)
p23=ggplot(scores_small_MAF_bin, aes(x=bin, y=ADS)) +
          geom_boxplot()+
          theme_minimal()+
          theme(text = element_text(family = 'Helvetica Neue'),
                plot.title = element_text(size = 19),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 15),
                panel.grid.major.x = element_blank() ,
                panel.grid.major.y = element_line( size=.1, color="grey85"))+
#          theme(axis.text.x  = element_text(angle=90, vjust=0.5))
          labs(title ="ADS distribution", subtitle = "MAF 0 to 0.05", x="Minor Allele Frequency (MAF)", y="Allele Dispersion Score (ADS)")

#Too big to save a svg
ggsave(filename = paste0("graph_output/",file_name,"_ADS_distrib_lowMAF.png"), p23, dpi =300 )

print("Number of variants per bin for supp data")
print(data_per_bin_small_MAF)

