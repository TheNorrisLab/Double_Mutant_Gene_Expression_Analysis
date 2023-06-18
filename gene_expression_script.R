#Script to take in gene data and output a sorted gene table with each gene and exp levels of double and single mutants

#This will output gene expression table.csv which will give all genes log2 change not sorted or ordered 
#gene exp dysregulation will give the ordered gene list of all the genes
#log2 change expression outside SE.csv will show the log2 change list only genes with expression is outside the standard error on the lfc

#Read csv files from computer 
#Comparing wild-type (N2) versus tdp-1;ceh-14 double mutant gene expression outputs from DESeq
setwd("C:/Users/my_directory/")
mut_1<-read.csv("N2_vs_tdp-1.csv", header=TRUE, row.names=1)
doublemut<-read.csv("N2vstdp1ceh14.csv",header=TRUE, row.names=1)
mut_2<-read.csv("N2vsceh14.csv",header=TRUE, row.names=1)

#names gives the names of mutant 1, mutant 2, doublemutant
names<-c('tdp-1', 'ceh-14', 'tdp-1;ceh-14')


#Make files with only WBGenes with double mutant padj values <.05 and not NA and sort the data
double_mut<-doublemut[!is.na(doublemut$padj)&(doublemut$padj<.05),]
doublemut_ordered<-double_mut[with(double_mut, order(double_mut[,2])), ]
WBGenes<-c(row.names(doublemut_ordered))
mut1<-mut_1[WBGenes,]
mut2<-mut_2[WBGenes,]

#Create a data table to hold the log2FoldChange for each mutant
gene_exp<-data.frame(mut1$log2FoldChange, mut2$log2FoldChange, doublemut_ordered$log2FoldChange)
colnames(gene_exp)<-names
row.names(gene_exp)<-WBGenes

#Construct the heatmap (needs to be a data matrix not a data frame)
datamat<-data.matrix(gene_exp)
row.names(datamat)<-WBGenes
heatmap( datamat, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))

#Save dataframe to a csv with significant gene expression values from double mutant
write.csv(gene_exp ,file="gene expression table.csv")

#Play with sorting the data in different ways and see how the new heatmap looks
datamatsort<-datamat[(order(datamat[,2])),]
heatmap( datamatsort, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(10,10))

#order gene exp data
#order from smallest difference either mutant has with the double mutant
#add a column witht the smaller difference between mutant and double mutant then order by that column
ordered_gene_exp<-cbind(gene_exp, 0, 0)
#give some genaric column names so we can access each column using dataframe$column command
colnames(ordered_gene_exp)<-c( 'mut1', 'mut2', 'dm', 'exp difference', 'exp not between single mutant exp')

#take out rows where either single mutants expression=NA
ordered_gene_exp<-ordered_gene_exp[!is.na(ordered_gene_exp$mut1) & !is.na(ordered_gene_exp$mut2),]

#define a function to check if a number is between two other numbers
isBetween = function(x, y, z){
  if ((x <= y && x >= z) || (x >=y && x <= z)){
    TRUE
  }
  else if (is.na(x) || is.na(y) || is.na(z)){
    TRUE
  }
  else{
    FALSE
  }
}


#Use a loop to find the minimum difference between single mutant and double mutant then check if it's between those values
for (i in 1:nrow(ordered_gene_exp)){ #for each row in ordered_gene_exp
  m<-min(abs(abs(ordered_gene_exp[i,1])-abs(ordered_gene_exp[i,3])), abs(abs(ordered_gene_exp[i,2])-abs(ordered_gene_exp[i,3])))
  ordered_gene_exp[i,4]<-m
  #if the expression is between the two single mutants, put false (0), else put true (1) in the last column
  if((ordered_gene_exp[i,3]> ordered_gene_exp[i,1] & ordered_gene_exp[i,3] < ordered_gene_exp[i,2]) | (ordered_gene_exp[i,3]> ordered_gene_exp[i,2] & ordered_gene_exp[i,3]<ordered_gene_exp[i,1])){
    ordered_gene_exp[i,5]<-FALSE
  }
  else{
    ordered_gene_exp[i,5]<-TRUE
  }
}
#Change the column names to match the specific genes you're using
colnames(ordered_gene_exp)<-c(names, 'exp difference', 'exp not between single mutant exp')

#order the data based on largest difference in gene expression
gene_exp_ordered<-ordered_gene_exp[with(ordered_gene_exp, order(-ordered_gene_exp[,4])), ]
#order so that true values are on top (values OUTSIE the gene expression of single mutants)
gene_exp_tf<-gene_exp_ordered[order(-gene_exp_ordered[,5]),]

#Output to a csv
write.csv(gene_exp_tf, file="Gene Exp Dysregulation.csv")

#Create a heatmap with all genes ordered and all genes with expression > 2
heatmap(data.matrix(gene_exp_tf[,c(1,2,3)]), Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10), main='All genes')
heatmap(data.matrix(gene_exp_tf[-which(gene_exp_tf$`exp difference`<2),c(1,2,3)]), Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10), main='Significant Genes')


#look at lfcSE and order by how far outside lfcsE the genes are: hopefully this will minimize large errors from 
gene_lfcSE <- data.frame(mut1$log2FoldChange, mut1$lfcSE, mut2$log2FoldChange, mut2$lfcSE, doublemut_ordered$log2FoldChange, doublemut_ordered$lfcSE, 0)
colnames(gene_lfcSE)<-c( "Mut1 lfc" , "mut1 lfcSE", "mut2 lfc", "mut2 lfcSE", "double mutant lfc", "double mutant lfcSE", "Outside lfcSE")
row.names(gene_lfcSE)<- row.names(mut1)
#get rid of NA by using only rows from previously filtered data
geneLfcSE <- gene_lfcSE[(row.names(gene_exp_tf)),]
geneLfcSE <- data.frame(geneLfcSE, gene_exp_tf$`exp difference`)



#for each row if gene exp 
for (i in 1:nrow(geneLfcSE)){
  #if it's in between the gene levels plus or minus expression
  if(isBetween(geneLfcSE[i,5]+ geneLfcSE[i,6], geneLfcSE[i,1] -geneLfcSE[i,2], geneLfcSE[i,3]+geneLfcSE[i,4]) || isBetween(geneLfcSE[i,5]- geneLfcSE[i,6], geneLfcSE[i,1] -geneLfcSE[i,2], geneLfcSE[i,3]+geneLfcSE[i,4]) || isBetween(geneLfcSE[i, 5] + geneLfcSE[i,6], geneLfcSE[i,1] +geneLfcSE[i,2], geneLfcSE[i,3]-geneLfcSE[i,4]) ||isBetween(geneLfcSE[i, 5] - geneLfcSE[i,6], geneLfcSE[i,1] +geneLfcSE[i,2], geneLfcSE[i,3]-geneLfcSE[i,4]) ){
    geneLfcSE[i,7]=FALSE
  }

  else {
    geneLfcSE[i,7]= TRUE
  }
}

geneLfcSE_outside <- geneLfcSE[geneLfcSE$Outside.lfcSE == 1,]
write.csv(geneLfcSE_outside, file = "log2 change expression outside SE.csv")


#Sort by p value


