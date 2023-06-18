#Script to specifically identify gene expression / splicing that is dysregulated specifically in double mutants

#Order gene expression and splicing from most different to least different
#decide a statistical cutoff
#return relevant splice sites and WBGenes 

#use A3SS_max, A5SS_max, cassette_df_new, and retainedintron_new generated from splicing script for alternative splicing data
#Upload file with gene expression log2foldchange
filename1<-'C:/Users/my_directory/datatable.csv'
gene_exp<-read.csv(file = filename1, header = TRUE)

#order gene exp data
#Ordering possibilities
#1. order from absolute value of difference in mbl-1 plus exc-7 if both are significant
  #set er= error of RNAseq data that we might reasonably expect two different numbers to be the same (start with 1)
er=1
  #If either single mutant's gene exp is within 'er' number of the double mutant, don't count it as being significantly disregulated (towards bottom of the list)
#2. order from smallest difference either mutant has with the double mutant
#Try option 2
  #add a column witht the smaller difference between mutant and double mutant then order by that column
ordered_gene_exp<-(gene_exp)
ordered_gene_exp<-cbind(ordered_gene_exp, 0, 0)
colnames(ordered_gene_exp)<-c('Gene', 'mbl-1', 'exc-7', 'mbl-1;exc-7', 'exp difference', 'exp not between single mutant exp')

for (i in 1:nrow(gene_exp)){
  if (is.na(gene_exp[i,2])){
    m<-(gene_exp[i,3]-gene_exp[i,4])
    ordered_gene_exp[i,2]<-0
  }
  else if(is.na(gene_exp[i,3])){
    m<-(gene_exp[i,2]-gene_exp[i,4])
    ordered_gene_exp[i,3]<-0
  }
  else {
    m<-min(abs(abs(gene_exp[i,2])-abs(gene_exp[i,4])), abs(abs(gene_exp[i,3])-abs(gene_exp[i,4])))
  }
  ordered_gene_exp[i,5]<-m
  if(abs(ordered_gene_exp[i,4])> abs(ordered_gene_exp[i,2]) & abs(ordered_gene_exp[i,4])> abs(ordered_gene_exp[i,3])){
    ordered_gene_exp[i,6]<-TRUE
  }
  else{
    ordered_gene_exp[i,6]<-FALSE
  }
}
gene_exp_ordered<-ordered_gene_exp[with(ordered_gene_exp, order(-ordered_gene_exp[,5])), ]
gene_exp_tf<-gene_exp_ordered[order(-gene_exp_ordered[,6]),]
gene_exp_cut<-gene_exp_tf[gene_exp_tf[,5]>er,]


#Order genes higher if they are not between the expression levels of either single mutant


#Splicing Significance A3SS_max

ordered_A3SS<-cbind(A3SS_max, 0, 0)
colnames(ordered_A3SS)<-c("Gene", "mbl-1;exc-7 PSI", "mbl-1 PSI", "exc-7 PSI", "PSI Difference", "PSI Outside either")

for (i in 1:nrow(ordered_A3SS)){
  ordered_A3SS[i,5]<-min(abs(ordered_A3SS[i,3]-ordered_A3SS[i,2]), abs(ordered_A3SS[i,4]-ordered_A3SS[i,2]))
  if((ordered_A3SS[i,2]<=ordered_A3SS[i,3] & ordered_A3SS[i,2]>=ordered_A3SS[i,4]) | (ordered_A3SS[i,2]>=ordered_A3SS[i,3] & ordered_A3SS[i,2]<=ordered_A3SS[i,4]))
  {ordered_A3SS[i,6]<-FALSE}
  else {
    ordered_A3SS[i,6]<-TRUE
  }
}

A3SS_order<-ordered_A3SS[with(ordered_A3SS, order(-ordered_A3SS$`PSI Difference`)),]
A3SS_order<-A3SS_order[order(-A3SS_order$`PSI Outside either`),]

#Splicing Significance A5SS_max

ordered_A5SS<-cbind(A5SS_max, 0, 0)
colnames(ordered_A5SS)<-c("Gene", "mbl-1;exc-7 PSI", "mbl-1 PSI", "exc-7 PSI", "PSI Difference", "PSI Outside either")

for (i in 1:nrow(ordered_A5SS)){
  ordered_A5SS[i,5]<-min(abs(ordered_A5SS[i,3]-ordered_A5SS[i,2]), abs(ordered_A5SS[i,4]-ordered_A5SS[i,2]))
  if((ordered_A5SS[i,2]<=ordered_A5SS[i,3] & ordered_A5SS[i,2]>=ordered_A5SS[i,4]) | (ordered_A5SS[i,2]>=ordered_A5SS[i,3] & ordered_A5SS[i,2]<=ordered_A5SS[i,4]))
  {ordered_A5SS[i,6]<-FALSE}
  else {
    ordered_A5SS[i,6]<-TRUE
  }
}

A5SS_order<-ordered_A5SS[with(ordered_A5SS, order(-ordered_A5SS$`PSI Difference`)),]
A5SS_order<-A5SS_order[order(-A5SS_order$`PSI Outside either`),]


#Splicing Significance cassette_df_new

ordered_cassette<-cbind(cassette_df_new, 0, 0)
colnames(ordered_cassette)<-c("Gene", "mbl-1;exc-7 PSI", "mbl-1 PSI", "exc-7 PSI", "PSI Difference", "PSI Outside either")

for (i in 1:nrow(ordered_cassette)){
  ordered_cassette[i,5]<-min(abs(ordered_cassette[i,3]-ordered_cassette[i,2]), abs(ordered_cassette[i,4]-ordered_cassette[i,2]))
  if((ordered_cassette[i,2]<=ordered_cassette[i,3] & ordered_cassette[i,2]>=ordered_cassette[i,4]) | (ordered_cassette[i,2]>=ordered_cassette[i,3] & ordered_cassette[i,2]<=ordered_cassette[i,4]))
  {ordered_cassette[i,6]<-FALSE}
  else {
    ordered_cassette[i,6]<-TRUE
  }
}

cassette_order<-ordered_cassette[with(ordered_cassette, order(-ordered_cassette$`PSI Difference`)),]
cassette_order<-cassette_order[order(-cassette_order$`PSI Outside either`),]

#Splicing Significance retainedintron_new

ordered_retainedintron<-cbind(retainedintron_new, 0, 0)
colnames(ordered_retainedintron)<-c("Gene", "mbl-1;exc-7 PSI", "mbl-1 PSI", "exc-7 PSI", "PSI Difference", "PSI Outside either")

for (i in 1:nrow(ordered_retainedintron)){
  ordered_retainedintron[i,5]<-min(abs(ordered_retainedintron[i,3]-ordered_retainedintron[i,2]), abs(ordered_retainedintron[i,4]-ordered_retainedintron[i,2]))
  if((ordered_retainedintron[i,2]<=ordered_retainedintron[i,3] & ordered_retainedintron[i,2]>=ordered_retainedintron[i,4]) | (ordered_retainedintron[i,2]>=ordered_retainedintron[i,3] & ordered_retainedintron[i,2]<=ordered_retainedintron[i,4]))
  {ordered_retainedintron[i,6]<-FALSE}
  else {
    ordered_retainedintron[i,6]<-TRUE
  }
}

retainedintron_order<-ordered_retainedintron[with(ordered_retainedintron, order(-ordered_retainedintron$`PSI Difference`)),]
retainedintron_order<-retainedintron_order[order(-retainedintron_order$`PSI Outside either`),]

#download the files
#write.csv(A3SS_order ,paste('C:/Users/my_directory/RNA seq splicing/','A3SS_ordered.csv',sep=""))
#write.csv(A5SS_order ,paste('C:/Users/my_directory/RNA seq splicing/','A5SS_ordered.csv',sep=""))
#write.csv(retainedintron_order ,paste('C:/Users/my_directory/RNA seq splicing/','retainedintron_ordered.csv',sep=""))
#write.csv(cassette_order ,paste('C:/Users/my_directory/RNA seq splicing/','cassette_ordered.csv',sep=""))