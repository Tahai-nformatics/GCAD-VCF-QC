#!/usr/bin/env Rscript

suppressMessages(
 library(dplyr)
)
library (optparse)

# parse arguments
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

# Read-in SNV files
summary_data = data.frame()

for (i in c(1:22)){

   f_name = paste("./chr", i, "/summary.snv.", subset ,".tsv", sep='')
   if(file_test('-f', f_name)){
   
     # read only specific cols (Mend_Incon VFLAGS VTYPE)
     stuff <- read.delim(f_name, colClasses = c(rep("NULL", 22), "numeric", rep("NULL",4), "character", rep("NULL", 5), "character"), header=T)
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
# Pre-QC, total # variants by type
col3_snv   <- df1 %>% filter(SNV==T) %>% count()
col3_indel <- df1 %>% filter(SNV==F) %>% count()

# Column 4
# Pre-QC, total # variants with mi-errors
col4_snv <- 0
col4_indel <- 0
if (df1 %>% filter(SNV==T, Mend_Incon > 0) %>% count() > 0){
    col4_snv   <- df1 %>% filter(SNV==T, Mend_Incon > 0) %>% select(Mend_Incon) %>% sum()
}
if (df1 %>% filter(SNV==F, Mend_Incon > 0) %>% count() > 0){
    col4_indel <- df1 %>% filter(SNV==F, Mend_Incon > 0) %>% select(Mend_Incon) %>% sum()
}

# Column 5
col5_snv <- sprintf("%0.4f%%", col4_snv / col3_snv * 100)
col5_indel <- sprintf("%0.4f%%", col4_indel / col3_indel * 100)

# Column 6
# Post-QC, total # variants (VFLAG==0)
col6_snv   <- df1 %>% filter(SNV==T, V0==T) %>% count()
col6_indel <- df1 %>% filter(SNV==F, V0==T) %>% count()

# Column 7
# Post-QC, total # variants (VFLAG==0) with mi-errors
col7_snv <- 0
col7_indel <- 0
if (df1 %>% filter(SNV==T, V0==T, Mend_Incon > 0) %>% count() > 0){
    col7_snv   <- df1 %>% filter(SNV==T, V0==T, Mend_Incon > 0) %>% select(Mend_Incon) %>% sum()
}
if (df1 %>% filter(SNV==F, V0==T, Mend_Incon > 0) %>% count() > 0){
    col7_indel <- df1 %>% filter(SNV==F, V0==T, Mend_Incon > 0) %>% select(Mend_Incon) %>% sum()
}

# Column 8
col8_snv <- sprintf("%0.4f%%", col7_snv / col6_snv * 100)
col8_indel <- sprintf("%0.4f%%", col7_indel / col6_indel * 100)


# Output Header
header <- paste("Dataset", "Variant Type",
                "N, Total Bi-Allelic Variants (pre-QC)",
                "N_MI, Total Variants with Error (pre-QC)",
                "% variants with Mendelian Error",
                "Total Passing Variants (post-QC)",
                "Total Passing with Mendelian error (post-QC)",
                "% passing variants with Mendelian error",
                sep="\t")

# Data output
output1 <- paste(subset, "SNV",  col3_snv,  col4_snv,   col5_snv,   col6_snv,   col7_snv,   col8_snv, sep="\t")
output2 <- paste(subset, "Indel",col3_indel,col4_indel, col5_indel, col6_indel, col7_indel, col8_indel, sep="\t")

cat(header, output1, output2, sep="\n")
