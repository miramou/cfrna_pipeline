#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

file_path=args[1]

ribo=read.table(file_path,header=FALSE)
ribo[3,1]=ribo[1,1]/ribo[2,1]
write.table(ribo, file_path, append=FALSE, col.names=FALSE, row.names=FALSE)
