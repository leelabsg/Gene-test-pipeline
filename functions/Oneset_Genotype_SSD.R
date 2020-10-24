#Oneset_genotype_SSD
Oneset_Genotype_SSD<-function(gene,anno,bfile,number){
  #gene you want, annovar file made by Read_Annovar, plink filename, initial number of genes in each setIDfile
  anno<-fread(anno,fill=T)
  genelist<-unique(anno$Gene.refGene)
  
  geneloc<-which(genelist==gene)
  print(geneloc)
  setnumber<-floor(geneloc/number)+1
  quotient<-floor(length(genelist)/number)
  if(setnumber!=quotient+1){
    File.SSD<-paste0('selected_',as.character(setnumber),bfile,'.SSD')
    File.Info<-paste0('selected_',as.character(setnumber),bfile,'.SSD.info')
  }else{
    File.SSD<-paste0('selected_last',bfile,'.SSD')
    File.Info<-paste0('selected_last',bfile,'.SSD.info')
  }
  SSD.INFO<-Open_SSD(File.SSD,File.Info)
  id<-SSD.INFO$SetInfo[which(SSD.INFO$SetInfo$SetID==gene),]$SetIndex
  Get_Genotypes_SSD(SSD.INFO,id)
}
