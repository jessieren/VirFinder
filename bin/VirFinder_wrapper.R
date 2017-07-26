#!/usr/bin/env Rscript

library("optparse")
 
option_list = list(
  make_option(c("-f", "--fasta"), type="character", default=NULL, 
              help="fasta file name", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="out.tsv", 
              help="output file name [default= %default]", metavar="character"),
	make_option(c("-s", "--sort") , action="store_true",default=FALSE,
              help="Sort by p-value (or q-value if enabled) [default= %default]"),
	make_option(c("-q", "--qvalue") , action="store_true",default=FALSE,
              help="Calculate q-values [default= %default]")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$fasta)){
  print_help(opt_parser)
  stop("Must specify input .fasta file", call.=FALSE)
}

library(VirFinder)

print(paste0("Running VirFinder on ", opt$fasta))

predResult <- VF.pred(opt$fasta)
if (opt$qvalue){
predResult$qvalue <- VF.qvalue(predResult$pvalue)
}
if (opt$sort && opt$qvalue){
   predResult<-predResult[order(predResult$qvalue),]
}else{
  predResult<-predResult[order(predResult$pvalue),]
}

write.table(file=opt$out, x=predResult, sep="\t", quote=FALSE, row.names=FALSE)

print(paste0("VirFinder was successfully run and saved to ", opt$out))