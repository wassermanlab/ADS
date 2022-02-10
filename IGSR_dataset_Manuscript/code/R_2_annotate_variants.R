##Created by S. Correard
#Updated Feb 9th, 2022

#This script is intended to do two graph
#1. The graph of the ADS per annotation
#2. The graphs of the ADS per MAF per annotation

#files
current_path_name=getwd()
input_file_name="1KG_dataset_27022019.GRCh38.phased"

#install.packages('gridExtra')

library(ggplot2)
library("ggplot2")
library(dplyr)
library(gridExtra)
library(svglite)

###################
scores=read.table(file=paste0("subset_table_score_",input_file_name, "_allchr_annotated.variant_function"), header = FALSE)
colnames(scores)=c("annotation", "gene", "chr", "pos1_1", "pos1_2", "ref", "alt", "pos_i", "V4", "V5", "n_ind_hom1_i", "n_ind_hom2_i", "n_ind_het_i", "n_ind_mis_i", "MAF_i", "ADS", "length_unique_pop_i", "length_unique_pop2_i")

head(scores)
nrow(scores)

#Annotation to keep for the graph
annot=c("exonic", "splicing", "UTR5", "UTR3", "intronic", "upstream", "downstream", "intergenic", "ncRNA_exonic",  "ncRNA_intronic", "upstream;downstream", "ncRNA_splicing")
#skipped_annot=c("ncRNA_exonic;splicing", "exonic;splicing", "UTR5;UTR3", "ncRNA_splicing")

###################

##One graph with all the annotations, no matter the MAF
p12=ggplot(scores, aes(x = reorder(annotation, ADS, FUN = median), y=ADS)) +
        geom_boxplot()+
	theme_minimal()+
	theme(text = element_text(family = 'Helvetica Neue'),
        plot.title = element_text(size = 18),
        axis.title = element_text(size = 20),
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line( size=.1, color="grey85"),
	axis.text.x  = element_text(angle=45, vjust=0.5))+
	labs(title ="ADS per annotation", y="Allele Dispersion Score (ADS)", x="Variants annotation")

ggsave(filename = paste0("graph_output/",input_file_name,"_ADS_distrib_by_annot.svg"), p12, dpi =200 )

#Print data in a table
data_per_annot <- count(scores, annotation)
data_per_annot_ordered=data_per_annot[order(data_per_annot$n),]
write.table(data_per_annot_ordered, file=paste0("graph_output/",input_file_name,"_data_per_annotation"), quote=FALSE, row.names = FALSE)

###################
#One graph for the annotations with more than 200 variants
#table with number of variant per annotation

annot_over_500=data_per_annot_ordered[which(data_per_annot_ordered$n>500),]
scores_annot_interest_table=subset(scores, annotation %in% annot_over_500$annotation)
write.table(scores_annot_interest_table, file=paste0("graph_output/",input_file_name,"_data_per_annotation_over500"), quote=FALSE, row.names = FALSE)


p13=ggplot(scores_annot_interest_table, aes(x = reorder(annotation, ADS, FUN = median), y=ADS)) +
	geom_boxplot()+
        theme_minimal()+
	        theme(text = element_text(family = 'Helvetica Neue'),
	        plot.title = element_text(size = 18),
        	axis.title = element_text(size = 20),
        	panel.grid.major.x = element_blank() ,
        	panel.grid.major.y = element_line( size=.1, color="grey85"),
        	axis.text.x  = element_text(angle=45, vjust=0.5))+
        scale_x_discrete(labels=c("splicing" = "Splicing", "exonic" = "Exonic", "UTR5" = "5' UTR", "UTR3" = "3' UTR", "ncRNA_splicing" = "Splicing of ncRNA", "upstream;downstream" = "Upstrean and/ndownstrean", "ncRNA_exnic" = "Exonic of ncRNA", "intronic" = "Intronic", "downstream" = "Downstream", "upstream" = "Upstream", "ncRNA_intronic" = "Intronic of ncRNA", "intergenic"="Intergenic"))+
        labs(title ="ADS per annotation", y="Allele Dispersion Score (ADS)", x="Variants annotation")
	
ggsave(filename = paste0("graph_output/",input_file_name,"_ADS_distrib_by_annot_over500.svg"), p13, dpi =200 )


###################
#One graph per annotation, with binned MAF
for (j in 1:length(annot)){
	print(j)
	annot_type=annot[j]
	print(annot_type)
	#Subset table for values in this MAF
	scores_annot_interest=which(scores$annotation==annot_type)
	scores_annot_interest_table=scores[scores_annot_interest,]
        scores_annot_interest_table_bin = scores_annot_interest_table %>% mutate( bin=cut_width(MAF_i, width=0.01, boundary=0) )

	#Calcule the number of observation per bin
	data_per_bin <- count(scores_annot_interest_table_bin, bin)
	write.table(data_per_bin, file=paste0("graph_output/",input_file_name,"_data_per_bin_",annot_type), quote=FALSE, row.names = FALSE)

	#Do the plot
	p3=ggplot(scores_annot_interest_table_bin, aes(x=bin, y=ADS)) +
		geom_boxplot()+
	  	theme_minimal()+
		theme(axis.text.x  = element_text(angle=45, vjust=0.5),
		        plot.title = element_text(size = 19),
        		axis.title = element_text(size = 19))+
		labs(title =paste0( "ADS distribution - ",annot_type), y="Allele Dispersion Score (ADS)", x="Minor Allele Frequency (MAF)")

	ggsave(filename = paste0("graph_output/",input_file_name,"_ADS_distrib_",annot_type,"_boxplot.svg"), p3, dpi =200 )

}


