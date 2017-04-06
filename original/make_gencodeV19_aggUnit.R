#############
##Goal
  #create a genic unit for ENSEMBLE basic set of genes such that each genic unit is a genomic region spanning all the transcripts of a specific gene
##output
  #  - note the output file is 1 based
  #  - the variable agg_start in output denotes the start (lower bound) of the genic unit and not necessarily the TSS
  #  - the variable agg_end in output denotes the end(upper bound) of the genic unit and not necessarily the TES
  #  - TES and TSS can be inferred based on the strand of the genic unit
#############
rm(list=objects())
library(GWASTools)
sessionInfo()
base.path<-"/projects/topmed"
gt<-read.table(file.path(base.path,"downloaded_data/Gencode/v19/gencode.v19.annotation.gtf.gz"),as.is=TRUE,sep="\t"); dim(gt)
colnames(gt)<-c("seqname","source","feature","start","end","score","strand","frame","attribute") #see http://uswest.ensembl.org/info/website/upload/gff.html#fields to figure the feilds of the gtf file

table(gt$feature,grepl("tag basic",gt$attribute))
table(gt$feature=="transcript",grepl("tag basic",gt$attribute))

# subset gtf to give a dataframe with transcripts for only the basic set of genes
gtb<-gt[gt$feature=="transcript" & grepl("tag basic",gt$attribute),]; dim(gtb)

#---------------
# extract values from attributes variables for easy manipualtion
#---------------
#get transcript_id for each transcript
gtb$transcript_id<-gsub("; gene_type.*","",(gsub(".*transcript_id.","", gtb$attribute)))
table(duplicated(gtb$transcript_id)) # these should be unique

#get transcript type
gtb$transcript_type<-gsub("; transcript_status.*","",(gsub(".*transcript_type.","", gtb$attribute)))

# get gene_name
gtb$gene_name<-gsub("; transcript_type.*","",(gsub(".*gene_name.","", gtb$attribute)))

#get gene_id for each transcript
gtb$gene_id<-gsub("; transcript_id.*","",(gsub(".*gene_id.","", gtb$attribute)))
table(duplicated(gtb$gene_id))

ugid<-unique(gtb$gene_id) ; length(ugid)

#---------------
# sanity check: check that a same gene_id is not present on diff chrs. 
#---------------
for (i in ugid){
ck<-gtb$seqname[gtb$gene_id==i]
ckk<-length(unique(ck))
stopifnot (ckk==1)
}

#---------------
# for a every unique gene_id, take the lowest start of all its transcript as agg_start and the highest end of all its transcript as its agg_end
#---------------
N<-length(ugid)
agf <- data.frame(chr=character(N),
                 strand=character(N),
                 gene_id=character(N),
                 gene_name=character(N),
                 agg_start=integer(N),
				 agg_end=integer(N),
				 source=character(N),
				 transcript_type=character(N),
				 merge.count=integer(N),
                 stringsAsFactors=FALSE)
				
agf$gene_id<-ugid

for (i in 1: nrow(agf)){
s<-gtb[gtb$gene_id==agf$gene_id[i],]
agf$gene_name[i]<-paste(unique(s$gene_name),collapse=",") 
stopifnot(!grep(",",agf$gene_name[i])) # expect this to be only one value
agf$chr[i]<-unique(s$seqname)
agf$strand[i]<-paste(unique(s$strand),collapse=",")
stopifnot(!grep(",",agf$strand[i])) # expect this to be only one value
agf$agg_start[i]<-min(s$start)
agf$agg_end[i]<-max(s$end)
agf$source[i]<-paste(unique(s$source),collapse=",")
agf$transcript_type[i]<-paste(unique(s$transcript_type),collapse=",")
agf$merge.count[i]<-nrow(s)
}
agf$agg_size<-agf$agg_end -agf$agg_start
agf<-agf[order(agf$chr,agf$agg_start,agf$gene_id),]
dim(agf) #[1] 44540    10
unlist(lapply(agf,function(x) sum(is.na(x)))) # all are zero
write.table(agf,file=file.path(base.path,"/gac_data/aggregation_units/gene_based/gencode_v19/gencode.v19.BasicGeneUnits.txt"),quote=FALSE,row.names=FALSE,sep="\t")

#---------------
# make a data dictionary
#---------------
names(agf)
# [1] "chr"             "strand"          "gene_id"         "gene_name"
# [5] "agg_start"       "agg_end"         "source"          "transcript_type"
# [9] "merge.count"     "agg_size"

vars<-names(agf)
type<-unlist(lapply(agf,class))
desc<-c("chromosome in genome build GRCh37/hg19", "DNA strand. values {+,-}","gene identifier with ENSG prefix","Gene name","lower bound of genic unit","upper bound of genic unit","comma seprated lists values of annotation sources of transcripts in a genic unit", "comma seprated values of biotypes of transcripts in a genic unit","number of transcripts merged in the genic unit","length of the genic unit derived as a difference of agg_end agg_start" )
ded<-data.frame("variables"=vars,"type"=type,"description"=desc,stringsAsFactors=FALSE)
write.table(ded,file=file.path(base.path,"/gac_data/aggregation_units/gene_based/gencode_v19/gencode.v19.BasicGeneUnits_DD.txt"),quote=FALSE,row.names=FALSE,sep="\t")


