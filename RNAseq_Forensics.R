setwd("/home/nikolay-oskolkov/TARGETWISE/Projects/HMGCR/RNAseq_muscle_liver/RNAseq_Forensics/GENE_SETS/")

tissues<-c("Adipose","Blood","Brain","Heart","Liver","Lung","Muscle","Pancreas","Skin")
gene_sets<-list()
for(i in tissues)
{
  gene_sets[[i]]<-readLines(paste0(i,".txt"))
}

#data<-read.delim("../HMGCR_two_samples_salmon_counts/Human_Muscle_JMS.txt",
#                  header=TRUE,check.names=FALSE,row.names=1,sep="\t")
#data<-read.delim("../HMGCR_two_samples_salmon_counts/Human_Blood.txt",
#                  header=TRUE,check.names=FALSE,row.names=1,sep="\t")
#data<-read.delim("../HMGCR_two_samples_salmon_counts/Human_Pancreatic_Islets.txt",
#                  header=TRUE,check.names=FALSE,row.names=1,sep="\t")
#data<-read.delim("../HMGCR_two_samples_salmon_counts/Human_Liver.txt",
#                  header=TRUE,check.names=FALSE,row.names=1,sep="\t")
#data<-read.delim("../HMGCR_two_samples_salmon_counts/Mouse_heart.txt",
#                  header=TRUE,check.names=FALSE,row.names=1,sep="\t")
#data<-read.delim("../HMGCR_two_samples_salmon_counts/Mouse_Liver_Muscle.txt",
#                  header=TRUE,check.names=FALSE,row.names=1,sep="\t")
#rownames(data)<-toupper(rownames(data))
#data<-read.delim("../HMGCR_two_samples_salmon_counts/Great_Apes_Muscle.txt",
#                  header=TRUE,check.names=FALSE,row.names=1,sep="\t")

data<-read.delim("/home/nikolay-oskolkov/TARGETWISE/Projects/HMGCR/RNAseq_muscle_liver/RNAseq_Forensics/github/merged_tissues_organisms_baseline.txt",header=TRUE,check.names=FALSE,row.names=1,sep="\t")

overlap_list<-list()
for(j in 1:ncol(data))
{
print(paste0("SAMPLE ",colnames(data)[j]))
top_expressed_genes_per_sample<-toupper(rownames(data)[order(data[,j],decreasing=TRUE)][1:100])
fraction_overlap <- sapply(gene_sets, function(gset) {
  length(intersect(gset, top_expressed_genes_per_sample)) / 
    length(top_expressed_genes_per_sample)})
overlap_list[[j]]<-fraction_overlap
}
overlap_df<-data.frame(Reduce(rbind,overlap_list))
rownames(overlap_df)<-colnames(data)
overlap_df<-overlap_df/rowSums(overlap_df)
overlap_df

library("RColorBrewer")
classes <- levels(factor(colnames(overlap_df)))
cols<-brewer.pal(n=length(unique(classes)),"Paired")
class_colors <- setNames(cols[seq_along(classes)], classes)
barplot(t(overlap_df),beside=FALSE,col=class_colors[colnames(overlap_df)],legend=TRUE,
        args.legend=list(x="topright",inset=.02,legend=colnames(overlap_df),
                         fill=class_colors[colnames(overlap_df)],cex=0.7),
        names.arg=rownames(overlap_df),ylab="Tissue Fraction",xlab="Samples Ground Truth",
        main="Predictions on new samples",cex.names=0.5,las=2)
