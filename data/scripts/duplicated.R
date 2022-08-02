#!/usr/bin/env Rscript

# clean environment, not really for script execution though, I guess?
options(warn=-1)
rm(list=ls())

# import libraries
library(optparse) |> suppressMessages()
library(rtracklayer) |> suppressMessages()

# take arguments
option_list = list(
  make_option(c("-f", "--gtf"), type="character", default=NULL, 
              help="input gtf file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.gtf", 
              help="output gtf file name [default= %default]", metavar="character"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Print extra output [default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# halt if no input
if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# import files
gtf = rtracklayer::import(opt$gtf)
gtf_df = as.data.frame(gtf)

# filter gtf to only include rows that are not the same in their start, end and type
filtered_gtf = gtf_df[!duplicated(gtf_df[,c("start", "end", "type")]), ]
#filtered_gtf = gtf_df

# How many rows were left out?
abs = dim(gtf_df)[1] - dim(filtered_gtf)[1]
rel = round(100 - (dim(filtered_gtf)[1] / dim(gtf_df)[1] * 100), 2)
output = sprintf("In the gtf that was supplied, %s entries were found. %s (%s %%) were found to be duplicated entries.",
        dim(gtf_df)[1], abs, rel)
if (opt$verbose == TRUE) {
  print(output)
}


# write the out file
rtracklayer::export(object = filtered_gtf, con = opt$out)
