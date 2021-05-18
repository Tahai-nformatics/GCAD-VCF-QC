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

# column 5
# Biallelic Variant Count
col5_snv   <- nrow(filter(summary_data, VTYPE == "SNV") )
col5_indel <- nrow(filter(summary_data, VTYPE != "SNV") )

# column 6
# Bi-allelic variants outside of target capture region: Total # of variants with VFLAG 11
col6_snv   <- df1 %>% filter(V11==T, SNV==T) %>% count()
col6_indel <- df1 %>% filter(V11==T, SNV==F) %>% count()

# column 7
# Biallelic Variants with VQSR Tranche>99.7%- Total # of variants with VFLAG 1 which do not have VFLAG 11 
col7_snv   <- df1 %>% filter(V1==T, V11==F, SNV==T) %>% count()
col7_indel <- df1 %>% filter(V1==T, V11==F, SNV==F) %>% count()

# column 8
# Variants with All Genotypes at DP<10 or GQ<20- Total # of variants with VFLAG 2 which do not have VFLAG 11 and/or 1 
col8_snv   <- df1 %>% filter(V2==T, V11==F, V1==F, SNV==T) %>% count()
col8_indel <- df1 %>% filter(V2==T, V11==F, V1==F, SNV==F) %>% count()

# column 9
# Failed Monomorphic filter- Total # of variants with VFLAG 3 which do not have VFLAG 11, 1, and 2 
col9_snv   <- df1 %>% filter(V3==T, V11+V1+V2==F, SNV==T) %>% count()
col9_indel <- df1 %>% filter(V3==T, V11+V1+V2==F, SNV==F) %>% count()

# column 10
# Failed Missing Rate filter- Total # of variants with VFLAG 4 which do not have VFLAG 11, 1, 2, 3
col10_snv   <- df1 %>% filter(V4==T, V11+V1+V2+V3==F, SNV==T) %>% count()
col10_indel <- df1 %>% filter(V4==T, V11+V1+V2+V3==F, SNV==F) %>% count()

# column 11
# Failed Mean Depth filter- Total # of variants with VFLAG 5 which do not have VFLAG 11, 1, 2, 3, 4
col11_snv   <- df1 %>% filter(V5==T, V11+V1+V2+V3+V4==F, SNV==T) %>% count()
col11_indel <- df1 %>% filter(V5==T, V11+V1+V2+V3+V4==F, SNV==F) %>% count()

# column 12
# After variant QC- Total # of variants with VFLAG0
col12_snv   <- df1 %>% filter(V0==T, SNV==T) %>% count()
col12_indel <- df1 %>% filter(V0==T, SNV==F) %>% count()


# Output Header
header <- paste("Dataset", "Variant Type", "Total Variants", "Multi-allelic Variants", "Biallelic Variants", 
                "Bi-allelic variants outside of target capture region",
                "Biallelic Variants with VQSR Tranche>99.7%",
                "Variants with All Genotypes at DP<10 or GQ<20",
                "Failed Monomorphic filter",
                "Failed Missing Rate filter",
                "Failed Mean Depth filter",
                "After variant QC",
                sep="\t")

# Data output
output1 <- paste(subset, "SNV",  col5_snv , 0, col5_snv, col6_snv, col7_snv, col8_snv, col9_snv, col10_snv, col11_snv, col12_snv, sep="\t")
output2 <- paste(subset, "Indel",col5_indel ,0,col5_indel, col6_indel, col7_indel, col8_indel, col9_indel, col10_indel, col11_indel, col12_indel, sep="\t")

cat(header, output1, output2, sep="\n")
