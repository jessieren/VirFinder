#!/usr/bin/env Rscript

library("optparse")
 
option_list = list(
  make_option(c("-f", "--fasta"), type="character", default=NULL, 
              help="fasta file name", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="out.tsv", 
              help="output file name [default= %default]", metavar="character")
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
write.table(file=opt$out, x=predResult, sep="\t")

print(paste0("VirFinder was successfully run and saved to ", opt$out))