##Var_Info function
Var_Info<-function(anno,result,order=T){
  anno<-read.csv(anno)
  #get list of genes and prepare empty dataset
  genelist<-unique(anno$Gene.refGene)
  df=data.frame(chr=NA,gene=NA,start=NA,end=NA)
  
  for(i in 1:length(genelist)){
    row<-c(0,0,0,0)
    row[1]<-anno[which(anno$Gene.refGene==genelist[i]),]$Chr[1]
    row[2]<-genelist[i]
    row[3]<-min(anno[which(anno$Gene.refGene==genelist[i]),]$Start)
    row[4]<-max(anno[which(anno$Gene.refGene==genelist[i]),]$Start)
    df<-rbind(df,row)
  }
  
  
  varinfo_table_a<-as.data.frame.matrix(table(anno$Gene.refGene,anno$Func.refGene))
  varinfo_table_b<-as.data.frame.matrix(table(anno$Gene.refGene,anno$ExonicFunc.refGene))
  varinfo_table<-cbind(varinfo_table_a,varinfo_table_b)
  varinfo_table<-tibble::rownames_to_column(varinfo_table,'gene')
  varinfo_table<-left_join(df,varinfo_table,by='gene')
  varinfo_table<-varinfo_table[2:nrow(varinfo_table),]
  
  varinfo_table$Lossoffunction<-varinfo_table$`frameshift deletion`+varinfo_table$`frameshift insertion`
  +varinfo_table$startloss+varinfo_table$stopgain+varinfo_table$stoploss
  
  write.table(varinfo_table, file=result,col.names=T,row.names=T)
  rm(varinfo_table)
  rm(varinfo_table_a)
  rm(varinfo_table_b) 
}
