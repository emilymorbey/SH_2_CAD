#Proxy matching

.libPaths("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/R_lib/")
sink("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/proxies_log_F_S.log")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

#Load required packages
install.packages("tidyverse", dependencies = TRUE)
library(tidyverse) 
install.packages("data.table")
library(data.table)

#### Read in the proxy file### 
Intermediate_Johns <- read.table(file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/M_Testosterone/merged_M_Testosterone_proxies_r06.txt", head=TRUE, sep="\t", row.names=NULL)
(head)
#Rename Proxy column into "SNP" (but keep proxy column as well)
Intermediate_Johns$SNP <- Intermediate_Johns$Proxy

#### Read in original exposure GWAS###
Input_data <- read.table(file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/GWAS_outputs/M_Testosterone_GWAS.txt", head=TRUE, sep="\t", row.names=NULL)

### Match the proxy list to the exposure GWAS SNPs (and their respective p-values)
Intermediate_SNPs <- merge(Intermediate_Johns, Input_data, by="SNP") 
head(Intermediate_SNPs)
#Only retain proxies that show a significant association with exposure - !!!!!!!might need to change p-value name!!!!!!! 
Intermediate_SNPs <- Intermediate_SNPs[which(Intermediate_SNPs$P_BOLT_LMM < 5e-5), ] 

print(Intermediate_SNPs)
#Convert the correlation value (R2) to numeric to be able to sort on this variable
Intermediate_SNPs$R2 <- as.numeric(as.character(Intermediate_SNPs$R2))

### we identify a list of variants that need proxies
look_up <- Intermediate_SNPs[!duplicated(Intermediate_SNPs$Target), ]  

### Read in outcome data 
Outcome_data <- read.table(file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/GWAS_outputs/CAD_GWAS.txt", head=T, sep="\t", row.names=NULL)

### Merge outcome data to proxy list by SNP or CHR&POS (BOLT output uses CHR and BP) - !!!!!!!!might need to change column names!!!!!!!!!!!
#Outcome_data$CHR <- Outcome_data$Chr
#Outcome_data$POS <- Outcome_data$Pos
#Outcome_data$SNP <- Outcome_data$MarkerName
search_out <- merge(Intermediate_SNPs, Outcome_data, by=c("CHR","BP")) 
head(search_out)
out_data <- c(colnames(search_out))

### We loop over each of the SNPs that might need a proxy, and for each we select the SNP with the best R2 value
for(i in seq(from = 1, to = length(look_up$Target), by = 1)){
  primary <- paste(look_up$Target[i])  
  subset <- search_out[which(search_out$Target==primary), ] 
  if (length(subset$SNP) == 0) {
    next
  }
  proxy <- subset[which(subset$R2 == max(subset$R2)), ] 
  print(proxy)
  x <- proxy[1, ]
  out_data <- rbind(out_data, x)
}

### We write out the data that matches our proxy data
write.table(out_data, file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/M_Testosterone_CAD_r06.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)