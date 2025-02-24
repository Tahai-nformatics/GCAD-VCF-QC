#!/usr/bin/env Rscript

suppressMessages(
 library(dplyr)
)
library (optparse)

# parset arguments
option_list <- list ( 
     make_option (c("-s","--subset"), help="subset name"),
     make_option (c("-d","--dir"), help="directory to scan")
     )

parser <-OptionParser(option_list=option_list)
arguments <- parse_args (parser, positional_arguments=TRUE)

opt <- arguments$options
args <- arguments$args

subset <- strsplit(opt$subset, ",")
wkdir <- unlist(strsplit(opt$dir, ","))
setwd(wkdir)

summary_data = data.frame()

for (i in c(1:22)){

   f_name = paste("./chr", i, "/summary.snv.", subset ,".tsv", sep='')
   if(file_test('-f', f_name)){
   
    # read only cols 6-10, 28,34
    stuff <- read.delim(f_name, colClasses = c(rep("NULL", 5), rep("character", 5), rep("NULL", 17), "character", rep("NULL", 5), "character"), header=T)
    summary_data <- rbind(summary_data, stuff)
    rm(stuff)
   }else{
     message(paste("missing", f_name))
   }
}

df1 <-  summary_data %>% 
        mutate(SNV=grepl("SNV", VTYPE), VTYPE=NULL) %>% 
        mutate(
               V0 = grepl("\\b0\\b",VFLAGS),
               V1 = grepl("\\b1\\b",VFLAGS),
               V2 = grepl("\\b2\\b",VFLAGS),
               V3 = grepl("\\b3\\b",VFLAGS),
               V4 = grepl("\\b4\\b",VFLAGS),
               V5 = grepl("\\b5\\b",VFLAGS),
               V11= grepl("\\b11\\b",VFLAGS),
               VFLAGS=NULL
               ) 

# Column 3
# N
col3_snv   <- df1 %>% filter(SNV==T) %>% count()
col3_indel <- df1 %>% filter(SNV==F) %>% count()

# remaining columns:
snv_cols <- df1 %>%  filter(SNV==T) %>% summarize(
   F00=sum(as.numeric(Fail00)),
   F01=sum(as.numeric(Fail01)),
   F11=sum(as.numeric(Fail11)),
   MIS=sum(as.numeric(Missing)),
   GFD=sum(as.numeric(GT_Failed)),
   V2S=sum(V2)
   )

ind_cols <- df1 %>%  filter(SNV==F) %>% summarize(
   F00=sum(as.numeric(Fail00)),
   F01=sum(as.numeric(Fail01)),
   F11=sum(as.numeric(Fail11)),
   MIS=sum(as.numeric(Missing)),
   GFD=sum(as.numeric(GT_Failed)),
   V2S=sum(V2)
   )

# Output Header
header <- paste("Dataset", "Variant Type",
                "N, Total Bi-Allelic Variants (pre-QC)",
                "0/0",
                "0/1",
                "1/1",
                "Missing",
                "Set to Missing",
                "N, Variants with all genotypes set to missing (post-genotype-level-QC)",
                sep="\t")

# Data output
output1 <- paste(subset, "SNV",  col3_snv,  snv_cols$F00, snv_cols$F01, snv_cols$F11, snv_cols$MIS, snv_cols$GFD, snv_cols$V2S, sep="\t")
output2 <- paste(subset, "Indel",col3_indel,ind_cols$F00, ind_cols$F01, ind_cols$F11, ind_cols$MIS, ind_cols$GFD, ind_cols$V2S, sep="\t")

cat(header, output1, output2, sep="\n")
