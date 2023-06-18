#RNAseq splicing

#upload files into dataframes with the rownames as the AS event ID 

#Path to folder where every file is
setwd("C:/Users/my_directory/fox-1;mbl-1 RNAseq")

#Comparing wild-type (N2) versus fox-1;mbl-1 double mutant gene expression outputs from JUM


#Type gene names in order: double mutant, mutant1, mutant2
names<-c('fox-1;mbl-1', 'fox-1', 'mbl-1')

#double mutant
dm_A3SS<-read.csv("mbl-1_fox-1_3SS.csv", header = TRUE, row.names=2, stringsAsFactors=FALSE)
dm_A5SS<-read.csv('mbl-1_fox-1_5SS.csv', header = TRUE, row.names=2, stringsAsFactors=FALSE)
dm_cassette<-read.csv('mbl-1_fox-1_cassette.csv', header = TRUE, row.names=2, stringsAsFactors=FALSE)
dm_retainedintron<-read.csv('mbl-1_fox-1_retainedintron.csv', sep=',', header = TRUE, row.names=2, stringsAsFactors = FALSE)
#Make sure that they are doubles not factors so we can look at them
dm_retainedintron$deltaPSI_control.treatment<-as.double((dm_retainedintron$deltaPSI_control.treatment))

#mutant1 fox-1 in this case
mut1_A3SS<-read.csv('fox-1_splicing_3SS.csv', header=TRUE, row.names=2, stringsAsFactors=FALSE)
mut1_A5SS<-read.csv('fox-1_splicing_5SS.csv', header=TRUE, row.names=2, stringsAsFactors=FALSE)
mut1_cassette<-read.csv('fox-1_splicing_cassette.csv',  header=TRUE, row.names=2, stringsAsFactors=FALSE)
mut1_retainedintron<-read.csv('fox-1_splicing_retainedintron.csv', header=TRUE, row.names=2, stringsAsFactors=FALSE)
mut1_retainedintron$deltaPSI_control.treatment<-as.double(mut1_retainedintron$deltaPSI_control.treatment)

#mutant 2 exc-7 in this case
mut2_A3SS<-read.table('AS_differential_JUM_output_A3SS_events_pvalue_1_final_simplified.txt', header=TRUE, row.names=2, stringsAsFactors=FALSE)
mut2_A5SS<-read.table('AS_differential_JUM_output_A5SS_events_pvalue_1_final_simplified.txt', header=TRUE, row.names=2, stringsAsFactors=FALSE)
mut2_cassette<-read.table('AS_differential_JUM_output_cassette_exon_events_pvalue_1_final_simplified.txt', header=TRUE, row.names=2, stringsAsFactors=FALSE)
mut2_retainedintron<-read.table('AS_differential_JUM_output_intron_retention_pvalue_1_final_simplified.txt', header=TRUE, row.names=2, stringsAsFactors=FALSE)
mut2_retainedintron$deltaPSI_control.treatment<-as.double(mut2_retainedintron$deltaPSI_control.treatment)

#Use only significant genes from double mutant, sort by p value in double mutant
dms_A3SS<-data.frame(dm_A3SS[-which(dm_A3SS$pvalue>.05),])
dms_A5SS<-data.frame(dm_A5SS[-which(dm_A5SS$pvalue>.05),])
dms_cassette<-data.frame(dm_cassette[-which(dm_cassette$pvalue>.05),])
dms_retainedintron<-data.frame(dm_retainedintron[-which(dm_retainedintron$pvalue>.05),])

#make dataframes with splicing event in double mutant and deltaPSI_control-treatment for each of 4 splicing types
A5SS_sites<-row.names(dms_A5SS)
A3SS_sites<-row.names(dms_A3SS)
cassette_sites<-row.names(dms_cassette)
retainedintron_sites<-row.names(dms_retainedintron)


A5SS_df<-data.frame(dms_A5SS$ï..Gene , dms_A5SS$deltaPSI_control.treatment, mut1_A5SS[A5SS_sites,]$deltaPSI_control.treatment, mut2_A5SS[A5SS_sites,]$deltaPSI_control.treatment)
row.names(A5SS_df)<-c(A5SS_sites)

A3SS_df<-data.frame(dms_A3SS$ï..Gene ,dms_A3SS$deltaPSI_control.treatment, mut1_A3SS[A3SS_sites,]$deltaPSI_control.treatment, mut2_A3SS[A3SS_sites,]$deltaPSI_control.treatment)
row.names(A3SS_df)<-c(A3SS_sites)

cassette_df<-data.frame(dms_cassette$ï..Gene,dms_cassette$deltaPSI_control.treatment, mut1_cassette[cassette_sites,]$deltaPSI_control.treatment, mut2_cassette[cassette_sites,]$deltaPSI_control.treatment)
row.names(cassette_df)<-c(cassette_sites)

retainedintron_df<-data.frame(dms_retainedintron$ï..Gene,dms_retainedintron$deltaPSI_control.treatment, mut1_retainedintron[retainedintron_sites,]$deltaPSI_control.treatment, mut2_retainedintron[retainedintron_sites,]$deltaPSI_control.treatment)
row.names(retainedintron_df)<-c(retainedintron_sites)

#get rid of NA values and Inf and INF
#Retained Intron
retainedintron_new<-retainedintron_df[!is.na(retainedintron_df$mut1_retainedintron.retainedintron_sites....deltaPSI_control.treatment)& !is.na(retainedintron_df$mut2_retainedintron.retainedintron_sites....deltaPSI_control.treatment),]
retainedintron_n<-retainedintron_new[!is.infinite(retainedintron_new$mut1_retainedintron.retainedintron_sites....deltaPSI_control.treatment) & !is.infinite(retainedintron_new$mut2_retainedintron.retainedintron_sites....deltaPSI_control.treatment),]
retainedintron_new<-retainedintron_n


#Parse dataframes so that they report the max value only for 5' and 3' Splice Site
#Get rid of NA and Inf values for single mutants
#cassette_df already has only one value in each point, so does retainedintron_df but retainedintron_df has null and inf values
cassette_df_new<-cassette_df[(!is.na(cassette_df$mut1_cassette.cassette_sites....deltaPSI_control.treatment)) & !is.na(cassette_df$mut2_cassette.cassette_sites....deltaPSI_control.treatment) & !is.infinite((cassette_df$mut1_cassette.cassette_sites....deltaPSI_control.treatment)) & !is.infinite(cassette_df$mut2_cassette.cassette_sites....deltaPSI_control.treatment),]

#convert to strings, parse into a list of strings, convert back to list of ints, take max

A5SS_max<-matrix(0, nrow=nrow(A5SS_df), ncol=ncol(A5SS_df))
A5SS_max<- data.frame(A5SS_max)
row.names(A5SS_max)<-row.names(A5SS_df)
colnames(A5SS_max)<-colnames(A5SS_df)
A5SS_max$dms_A5SS.ï..Gene<-A5SS_df$dms_A5SS.ï..Gene
for (i in 1:nrow(A5SS_df)){
  for (j in 2:ncol(A5SS_df)){
    list1<-A5SS_df[i,j]
    lis<-toString(list1)
    spl<-strsplit(lis, ";")
    fin<-type.convert(spl, na.strings="NA", as.is=FALSE, dec=".", numerals="allow.loss")
    n<-as.numeric(unlist(fin))
    val<-max(n)
    if (val!=Inf & !is.na(val)){
      A5SS_max[i,j]=val
    }
    else{
      A5SS_max[i,j]=A5SS_max[i,j]
    }
  }
}

A3SS_max<-matrix(0, nrow=nrow(A3SS_df), ncol=ncol(A3SS_df))
A3SS_max<- data.frame(A3SS_max)
row.names(A3SS_max)<- row.names(A3SS_df)
colnames(A3SS_max)<-colnames(A3SS_df)
A3SS_max$dms_A3SS.ï..Gene<-A3SS_df$dms_A3SS.ï..Gene
for (i in 1:nrow(A3SS_df)){
  for (j in 2:ncol(A3SS_df)){
    list1<-A3SS_df[i,j]
    lis<-toString(list1)
    spl<-strsplit(lis, ";")
    fin<-type.convert(spl, na.strings="NA", as.is=FALSE, dec=".", numerals="allow.loss")
    n<-as.numeric(unlist(fin))
    val<-max(n)
    if (val!=Inf & !is.na(val)){
      A3SS_max[i,j]=val
    }
    else{
      A3SS_max[i,j]=A3SS_max[i,j]
    }
  }
}
A3SS_max<-A3SS_max[(order(A3SS_max[,1])),]

#Take the minimum difference between single mutant and double mutant expression and order them

#export data to csv

write.csv(retainedintron_new, file='RetainedIntron.csv')
write.csv(A5SS_max ,file= 'A5SS.csv')
write.csv(cassette_df ,file='cassette.csv')
write.csv(A3SS_max ,file='A3SS.csv')

#Splicing Significance 

#A3SS_max

ordered_A3SS<-cbind(A3SS_max, 0, 0)
colnames(ordered_A3SS)<-c("Gene", names, "PSI Difference", "PSI Outside either")

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
colnames(ordered_A5SS)<-c("Gene", names, "PSI Difference", "PSI Outside either")

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
#get rid of na values in the single mutant rows
ordered_cassette<-ordered_cassette[!is.na(ordered_cassette$mut1_cassette.cassette_sites....deltaPSI_control.treatment),]
ordered_cassette<-ordered_cassette[!is.na(ordered_cassette$mut2_cassette.cassette_sites....deltaPSI_control.treatment),]
colnames(ordered_cassette)<-c("Gene", names, "PSI Difference", "PSI Outside either")


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
ordered_retainedintron<-data.frame(retainedintron_new, 0, 0)
colnames(ordered_retainedintron)<-c("Gene", names, "PSI Difference", "PSI Outside either")

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

#write csv files
write.csv(A3SS_order ,file='A3SS_ordered.csv')
write.csv(A5SS_order, file='A5SS_ordered.csv')
write.csv(retainedintron_order ,file='retainedintron_ordered.csv')
write.csv(cassette_order, file='cassette_ordered.csv')