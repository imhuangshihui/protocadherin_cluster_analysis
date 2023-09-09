#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(data.table)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

PCD45 <- readRDS("/mnt/hwstor9k_data1/junehuang/PCD/beta.quan.rmSexProbes.rds")

# Create a list of 45 PCD genes
gene <- read.delim("/mnt/hwstor9k_data1/junehuang/PCD/PCD45GeneList.txt", header = FALSE, sep = "\n")
genelist = c(gene)

# Search gene listed in anno850k dataset
# And extract the probes of 45 PCD genes in anno850k dataset
probes <- rownames(ann850k[ann850k$UCSC_RefGene_Name %in% genelist$V1,])

# Extract the probe values of 45 PCD genes of 0/16 FYU from beta file
allrelatedprobes <- PCD45[rownames(PCD45) %in% probes, ]
df <- data.frame(allrelatedprobes)

# T-test between 0 YFU and 16 FYU for each probe
result=c()
for (i in rownames(allrelatedprobes)){
  FYU0 <- df[i, seq(0,94,by=2)]
  FYU16 <- df[i, seq(1,94,by=2)]
  Pvalue <- t.test(FYU0, FYU16, alternative = c("two.sided"))$p.value
  result= rbind(result, cbind(i, data.frame(Pvalue), ann850k[i, 'UCSC_RefGene_Name'])) 
}
#head(result, n=10)

write.table(result, file = "pvalueof45pcdgenes.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)

# T-test between 0 YFU and 16 YFU for each patient
patient=c()
for (i in seq(1, 93, by=2)){
  p <- t.test(allrelatedprobes[,i], allrelatedprobes[,i+1])$p.value
  patient=rbind(patient, cbind(colnames(allrelatedprobes)[i], data.frame(p)))
}
patient
write.table(patient, file = "PvalForEachPatient.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)


