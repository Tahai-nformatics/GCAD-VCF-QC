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
   
    # read only cols 28,34
    stuff <- read.delim(f_name, colClasses = c(rep("NULL", 27), "character",rep("NULL", 5), "character"), header=T)
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

# column 3
# All variants outside of Target Capture Region, vlflag=11
col3_snv   <- df1 %>% filter(V11==T, SNV==T) %>% count()
col3_indel <- df1 %>% filter(V11==T, SNV==F) %>% count()

# column 4
# All variants with VQSR Tranche>99.7, vlflag=1
col4_snv   <- df1 %>% filter(V1==T, SNV==T) %>% count()
col4_indel <- df1 %>% filter(V1==T, SNV==F) %>% count()

# column 5
# All variants with Genotypes at DP<10 | GQ<20, vlflag=2
col5_snv   <- df1 %>% filter(V2==T, SNV==T) %>% count()
col5_indel <- df1 %>% filter(V2==T, SNV==F) %>% count()

# column 6
# All variants failing Monomorphic filter, vlflag=3
col6_snv   <- df1 %>% filter(V3==T, SNV==T) %>% count()
col6_indel <- df1 %>% filter(V3==T, SNV==F) %>% count()

# column 7
# All variants failing Call Rate filter, vlflag=4
col7_snv   <- df1 %>% filter(V4==T, SNV==T) %>% count()
col7_indel <- df1 %>% filter(V4==T, SNV==F) %>% count()

# column 8
# All variants failing Mean Depth filter, vlflag=5
col8_snv   <- df1 %>% filter(V5==T, SNV==T) %>% count()
col8_indel <- df1 %>% filter(V5==T, SNV==F) %>% count()

# Output Header
header <- paste("Dataset", "Variant Type", 
                "All variants outside of Target Capture Region",
                "All variants with VQSR Tranche>99.7",
                "All variants with Genotypes at DP<10 or GQ<20",
                "All variants failing Monomorphic filter",
                "All variants failing Call Rate filter",
                "All variants failing Mean Depth filter",
                sep="\t")

# Data output
output1 <- paste(subset, "SNV",  col3_snv , col4_snv, col5_snv, col6_snv, col7_snv, col8_snv, sep="\t")
output2 <- paste(subset, "Indel",col3_indel,col4_indel,col5_indel,col6_indel,col7_indel,col8_indel, sep="\t")

cat(header, output1, output2, sep="\n")
