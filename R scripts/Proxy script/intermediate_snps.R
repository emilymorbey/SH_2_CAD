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
Intermediate_Johns <- read.table(file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/F_SHBG/F_SHBG_proxies_r06.txt", head=TRUE, sep="\t", row.names=NULL)

#Rename Proxy column into "SNP" (but keep proxy column as well)
Intermediate_Johns$SNP <- Intermediate_Johns$Proxy

#### Read in original exposure GWAS###
Input_data <- read.table(file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/GWAS_outputs/F_SHBG_GWAS.txt", head=TRUE, sep="\t", row.names=NULL)

### Match the proxy list to the exposure GWAS SNPs (and their respective p-values)
Intermediate_SNPs <- merge(Intermediate_Johns, Input_data, by="SNP") 
head(Intermediate_SNPs)
#Only retain proxies that show a significant association with exposure - !!!!!!!might need to change p-value name!!!!!!! 
Intermediate_SNPs_sig <- Intermediate_SNPs[which(Intermediate_SNPs$P_BOLT_LMM < 5e-5), ] 
Intermediate_SNPs_insig <- Intermediate_SNPs[which(Intermediate_SNPs$P_BOLT_LMM > 5e-5), ]



### We write out the data that matches our proxy data
write.table(Intermediate_SNPs_insig, file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Proxies/intermediate_SNPs_insig.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)