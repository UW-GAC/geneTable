rm(list=objects())
library(GWASTools)
sessionInfo()
library(plyr)
base.path<-"/projects/topmed"
gt<-read.table(file.path(base.path,"downloaded_data/Gencode/v19/gencode.v19.annotation.gtf.gz"),as.is=TRUE,sep="\t"); dim(gt)
colnames(gt)<-c("seqname","source","feature","start","end","score","strand","frame","attribute") #see http://uswest.ensembl.org/info/website/upload/gff.html#fields to figure the feilds of the gtf file

#---------------
# extract values from attributes variables for easy manipualtion
#---------------
#get transcript_id for each transcript
gt$transcript_id<-gsub("; gene_type.*","",(gsub(".*transcript_id.","", gt$attribute)))
table(duplicated(gt$transcript_id)) # these should be unique

#get transcript type
gt$transcript_type<-gsub("; transcript_status.*","",(gsub(".*transcript_type.","", gt$attribute)))
                 1430                              11202
# get gene_name
gt$gene_name<-gsub("; transcript_type.*","",(gsub(".*gene_name.","", gt$attribute)))

#get gene_id for each transcript
gt$gene_id<-gsub("; transcript_id.*","",(gsub(".*gene_id.","", gt$attribute)))


#get gene type type
gt$gene_type<-gsub("; gene_status.*","",(gsub(".*gene_type.","", gt$attribute)))



#####################
# get genes in basic set
#####################
gen<-gt[gt$feature=="gene",]; dim(gen) #57820    14 
gen.out<-gen[,c("seqname", "start", "end","strand","gene_id" , "gene_name","gene_type")]
names(gen.out)[names(gen.out)=="seqname"]<-"chr"

write.table(gen.out,file=file.path(base.path,"/gac_data/aggregation_units/gene_based/gencode_v19/gencode_v19_genes.txt"),quote=FALSE,row.names=FALSE,sep="\t")

#---------------
# make a data dictionary
#---------------
names(gen.out)
[1] "chr"       "start"     "end"       "strand"    "gene_id"   "gene_name"
[7] "gene_type"

vars<-names(gen.out)
type<-unlist(lapply(gen.out,class))
desc<-c("chromosome in genome build GRCh37/hg19","5' start position of the gene in genome build GRCh37/hg19","3' end position of the gene in genome build GRCh37/hg19", "DNA strand. values {+,-}","gene identifier with ENSG prefix","Gene name"," biotype of the gene")
ded<-data.frame("variables"=vars,"type"=type,"description"=desc,stringsAsFactors=FALSE)
write.table(ded,file=file.path(base.path,"/gac_data/aggregation_units/gene_based/gencode_v19/gencode_v19_genes_DD.txt"),quote=FALSE,row.names=FALSE,sep="\t")