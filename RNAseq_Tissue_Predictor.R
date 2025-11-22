setwd("/home/nikolay/WABI/O_Hansson/Gtex/")

############################################ DATA PREPROCESSING ###############################################
library("RColorBrewer")
annot<-read.delim("GTEx_v7_Annotations_SampleAttributesDS.txt",header=TRUE,comment.char="#",sep="\t")
tissues_to_select<-c("Blood","Brain","Skin","Adipose Tissue","Muscle","Heart","Lung","Pancreas","Liver")
annot<-annot[grepl(paste0("^(", paste(tissues_to_select, collapse="|"), ")$"),as.character(annot$SMTS)),]
annot<-annot[!grepl("Cells",as.character(annot$SMTSD)),]
cols<-brewer.pal(n=length(unique(annot$SMTS)),"Paired")

library("data.table")
expr<-as.data.frame(fread("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz",header=TRUE))
rownames(expr)<-paste0(expr$Name,"_",expr$Description)
expr$Name<-NULL; expr$Description<-NULL

intersect_samples<-intersect(as.character(annot$SAMPID),colnames(expr))
annot<-annot[match(intersect_samples,as.character(annot$SAMPID)),]
expr<-subset(expr,select=intersect_samples)
expr<-expr[rowMeans(expr)>3000,]
rownames(expr)<-toupper(gsub("_","",substr(rownames(expr),19,nchar(rownames(expr)))))

my_newdata<-read.delim("HMGCR_two_samples_salmon_counts/merged_tissues_organisms.txt",
                       header=TRUE,check.names=FALSE,sep="\t")
my_newdata<-as.data.frame(t(my_newdata))
my_newdata[1:2,1:10]

expr<-expr[rownames(expr)%in%colnames(my_newdata),]

test_annot<-do.call(rbind,lapply(tissues_to_select, function(tissue){head(annot[annot$SMTS==tissue,],5)}))
test_expr<-subset(expr,select=as.character(test_annot$SAMPID))

expr<-expr[,!colnames(expr)%in%as.character(test_annot$SAMPID)]
annot<-annot[!as.character(annot$SAMPID)%in%as.character(test_annot$SAMPID),]


######################################## DIMENSIONALITY REDUCTION ############################################
PC <- prcomp(t(log10(expr + 1)), center=TRUE, scale=FALSE)
expl_var <- PC$sdev^2/sum(PC$sdev^2)
barplot(expl_var[1:20],ylab="EXPLAINED VARIANCE",main="VARIANCE EXPLAINED BY PRINCIPAL COMPONENTS",
        names.arg=paste0("PC",seq(1:20)),col="darkgreen")
plot(PC$x[,1:2],main="GTEX HUMAN TISSUES PCA PLOT",xlab="PC1",ylab="PC2",col=cols[as.factor(annot$SMTS)],pch=19,cex=0.8)
legend("bottomleft", inset=.02, levels(as.factor(annot$SMTS)), fill=cols[as.factor(levels(as.factor(annot$SMTS)))])

library("Rtsne")
optPerp<-round(sqrt(dim(expr)[2]),0)
tsne_opt_perp <- Rtsne(t(log10(expr + 1)),initial_dims=20,verbose=TRUE,check_duplicates=FALSE,
                       perplexity=optPerp,dims=2,max_iter=1000,partial_pca=TRUE)
plot(tsne_opt_perp$Y,main="GTEX HUMAN TISSUES TSNE PLOT",xlab="tSNE1",ylab="tSNE2",col=cols[as.factor(annot$SMTS)],
     pch=19,cex=0.8)
legend("topleft", inset=.02, levels(as.factor(annot$SMTS)), fill=cols[as.factor(levels(as.factor(annot$SMTS)))])


######################################### TRAINIG MACHINE LEARNING #############################################
library("rpart"); library("rpart.plot")
X<-as.data.frame(t(expr)); y<-annot$SMTS
df <- data.frame(y = y, X)
fit <- rpart(y ~ ., data = df, method = "class")
rpart.plot(fit)

newdata<-as.data.frame(t(test_expr))
colnames(newdata)<-make.names(colnames(newdata))
probs<-predict(fit,newdata=newdata,type="prob")
print(probs)

f <- as.factor(test_annot$SMTS)
classes <- levels(f); tsne_colors <- cols[as.integer(f)]
class_colors <- setNames(cols[seq_along(classes)], classes)
barplot(t(probs),beside=FALSE,col=class_colors[colnames(probs)],legend=TRUE,
        args.legend=list(x="topright",inset=.02,legend=colnames(probs),fill=class_colors[colnames(probs)],cex=0.7),
        names.arg=test_annot$SMTS,ylab="Tissue Fraction",xlab="Samples Ground Truth",
        main="Predictions on reserved GTEX samples",cex.names=0.7,las=2)


################################# PREDICTION ON COMPLETELY DIFFERENT SAMPLES #################################
train_genes <- colnames(X)
my_newdata <- my_newdata[, train_genes]
colnames(my_newdata) <- make.names(colnames(my_newdata))

probs<-predict(fit,newdata=my_newdata,type="prob")
print(probs)
barplot(t(probs),beside=FALSE,col=class_colors[colnames(probs)],legend=TRUE,
        args.legend=list(x="topright",inset=.02,legend=colnames(probs),fill=class_colors[colnames(probs)],cex=0.7),
        names.arg=rownames(probs),ylab="Tissue Fraction",xlab="Samples Ground Truth",
        main="Predictions on new samples",cex.names=0.4,las=2)

write.table(probs,file="HMGCR_two_samples_salmon_counts/probs_tissues_organisms.txt",
            col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
