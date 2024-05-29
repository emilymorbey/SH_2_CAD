# MATCHING CHR:POS TO RSID FOR SNPS THAT EXIST IN SHBG AND DBP
# F_SHBG
library(readxl)
library(tidyverse)

# Specify your custom headers
custom_headers1 <- c("chrpos")
custom_headers2 <- c("rsid", "chr", "bp", "chrpos")

# Read the Excel file and use custom headers
chrpos_dat <- read_excel("RSID_CHRPOS_MATCHING/F_SHBG_SNPs_in_CAD_chrpos.xlsx", col_names = custom_headers)
rsid_dat <- read_excel("RSID_CHRPOS_MATCHING/F_SHBG_RSIDs_from_SHBG_GWAS.xlsx", col_names = custom_headers2)

combined <- merge(chrpos_dat, rsid_dat, by = 'chrpos')


library(openxlsx)

# Specify the file path where you want to save the Excel file
excel_file_path <- "RSID_CHRPOS_MATCHING/"

# Write the merged data frame to an Excel file
write.xlsx(combined, excel_file_path)










# M SHBG

library(readxl)
library(tidyverse)

# Specify your custom headers
custom_headers1 <- c("chrpos")
custom_headers2 <- c("rsid", "chr", "bp", "chrpos")

# Read the Excel file and use custom headers
chrpos_dat <- read_excel("RSID_CHRPOS_MATCHING/M_SHBG_SNPs_in_CAD_chrpos.xlsx", col_names = custom_headers)
rsid_dat <- read_excel("RSID_CHRPOS_MATCHING/M_SHBG_RSIDs_from_SHBG_GWAS.xlsx", col_names = custom_headers2)

combined <- merge(chrpos_dat, rsid_dat, by = 'chrpos')


library(openxlsx)

# Specify the file path where you want to save the Excel file
excel_file_path <- "RSID_CHRPOS_MATCHING/"

# Write the merged data frame to an Excel file
write.xlsx(combined, excel_file_path)





# F Testosterone

library(readxl)
library(tidyverse)

# Specify your custom headers
custom_headers1 <- c("chrpos")
custom_headers2 <- c("rsid", "chr", "bp", "chrpos")

# Read the Excel file and use custom headers
chrpos_dat <- read_excel("RSID_CHRPOS_MATCHING/F_T_SNPs_in_CAD_chrpos.xlsx", col_names = custom_headers)
rsid_dat <- read_excel("RSID_CHRPOS_MATCHING/F_T_RSIDs_from_GWAS.xlsx", col_names = custom_headers2)

combined <- merge(chrpos_dat, rsid_dat, by = 'chrpos')


library(openxlsx)

# Specify the file path where you want to save the Excel file
excel_file_path <- "RSID_CHRPOS_MATCHING/"

# Write the merged data frame to an Excel file
write.xlsx(combined, excel_file_path)








# m testosterone

library(readxl)
library(tidyverse)

# Specify your custom headers
custom_headers1 <- c("chrpos")
custom_headers2 <- c("rsid", "chr", "bp", "chrpos")

# Read the Excel file and use custom headers
chrpos_dat <- read_excel("RSID_CHRPOS_MATCHING/M_T_SNPs_in_CAD_chrpos.xlsx", col_names = custom_headers)
rsid_dat <- read_excel("RSID_CHRPOS_MATCHING/M_T_RSIDs_from_T_GWAS.xlsx", col_names = custom_headers2)

combined <- merge(chrpos_dat, rsid_dat, by = 'chrpos')


library(openxlsx)

# Specify the file path where you want to save the Excel file
excel_file_path <- "RSID_CHRPOS_MATCHING/"

# Write the merged data frame to an Excel file
write.xlsx(combined, excel_file_path)


