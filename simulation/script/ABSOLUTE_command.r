#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
input_loc <- args[1]
out_loc <- args[2]
input_name <- args[3]
cnv_input <- paste(input_loc,input_name,".addNeutral.cnv.txt",sep="")
snv_input <- paste(input_loc,input_name,".maf.txt",sep="")
out <- paste(out_loc,input_name,".R.txt",sep="")
library(ABSOLUTE)
RunAbsolute(cnv_input,0, 0.015,0.95, 10, NA,'SNP_6.0', input_name, out_loc, 1500, 0.05, 0.005, 'total', snv_input, 0)
load(paste(out_loc,input_name,".ABSOLUTE.RData",sep=""))
alpha <- seg.dat[["mode.res"]][["mode.tab"]][,"alpha"]
tau <- round(seg.dat[["mode.res"]][["mode.tab"]][,"tau"])

sink(out)
cat(substr(input_name,10,nchar(input_name)))
cat("\t")
cat(round(alpha[tau==2][1],2))
cat("\n")
sink()

