library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(data.table)
library(tidyverse)
library(tibble)
library(ggsignif)
library(ggpubr)

ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
beta <- readRDS("/mnt/hwstor9k_data1/junehuang/PCD/beta.quan.rmSexProbes.rds")

foxj1 <- data.frame(beta[rownames(beta)] %in% )


# 從UCSC browser中挑出來1156個CpGSite, 根據這些site從beta文件中找對應的病人的site只有876個
pcdh_site = read.delim("pcdh.site.bed", header = FALSE, sep = "\t")
# Find all the sites for each clusters no matter
for (i in pcdh_site$V4){
  aa=rbind(aa, cbind(pcdh_site[pcdh_site$V4 == i, 1:3],i, 
                     ann850k[i, 'UCSC_RefGene_Name']))
}
value <- data.frame(beta[rownames(beta) %in% pcdh_site$V4, ])

# Create the vectors for column names
conname_0YFU <- scan("/home/junehuang/pcd/control_0YFU.txt", what = character(0))
conname_16YFU <- scan("/home/junehuang/pcd/control_16YFU.txt", what = character(0))
decname_0YFU <- scan("/home/junehuang/pcd/decline_0YFU.txt", what = character(0))
decname_16YFU <- scan("/home/junehuang/pcd/decline_16YFU.txt", what = character(0))

# Create the sub-dataset for each catagory
contr_0YFU <- value[, conname_0YFU]
contr_16YFU <- value[, conname_16YFU]
dec_0YFU <- value[, decname_0YFU]
dec_16YFU <- value[, decname_16YFU]

# Combine the multiple col values into one col and create the data frame
c0 <- c(contr_0YFU[,1])
for(i in 2:23){
  c0 <- c(c0,contr_0YFU[,i])
}
c16 <- c(contr_16YFU[,1])
for(i in 2:23){
  c16 <- c(c16, contr_16YFU[,i])
}
d0 <- c(dec_0YFU[,1])
for(i in 2:24){
  d0 <- c(d0, dec_0YFU[,i])
}
d16 <- c(dec_16YFU[,1])
for(i in 2:24){
  d16 <- c(d16, dec_16YFU[,i])
}
###################################DATASET######################################
my_data <- data.frame(
  val = c(c0, c16, d0, d16),
  diagnosis = rep(c("control(n=23)", "decline(n=24)"), times = c(40296, 42048)),
  year = rep(c("0YFU", "16YFU", "0YFU", "16YFU"), times = c(20148, 20148, 21024, 21024))
)
control <- my_data[my_data$diagnosis == "control(n=23)", ]
decline <- my_data[my_data$diagnosis == "decline(n=24)", ]
base <- my_data[my_data$year == "0YFU", ]
later <- my_data[my_data$year == "16YFU", ]

############################STATISTICS#########################################
group_by(control, year) %>%
  summarise(
    count = n(),
    mean = mean(val, na.rm = TRUE),
    sd = sd(val, na.rm = TRUE)
  )
t.test(val ~ year, data = control,
       alternative = "two.sided", paired = TRUE)

group_by(decline, year) %>%
  summarise(
    count = n(),
    mean = mean(val, na.rm = TRUE),
    sd = sd(val, na.rm = TRUE)
  )
t.test(val ~ year, data = decline,
       alternative = "two.sided", paired = TRUE)

group_by(base, diagnosis) %>%
  summarise(
    count = n(),
    mean = mean(val, na.rm = TRUE),
    sd = sd(val, na.rm = TRUE)
  )
t.test(val ~ diagnosis, data = base, alternative = "two.sided")

group_by(later, diagnosis) %>%
  summarise(
    count = n(),
    mean = mean(val, na.rm = TRUE),
    sd = sd(val, na.rm = TRUE)
  )
t.test(val ~ diagnosis, data = later, alternative = "two.sided")
############################PLOT###############################################
p1 <- ggboxplot(my_data, x = "year", y = "val",
                color = "diagnosis", palette = "jco")
p1 + stat_compare_means(aes(group = diagnosis), method = "t.test")

p2 <- ggboxplot(my_data, x = "diagnosis", y = "val",
                color = "year", palette = "jco")
p2 + stat_compare_means(aes(group = year), method = "t.test", paired = T)

##########################FOR EACH SIET#########################################
aa=c()
for (i in rownames(value)){
  control0 <- value[i, conname_0YFU]
  control16 <- value[i, conname_16YFU]
  decline0 <- value[i, decname_0YFU]
  decline16 <- value[i, decname_16YFU]
  p_0YFU <- t.test(control0, decline0, alternative = c("two.sided"))$p.value
  p_16YFU <- t.test(control16, decline16, alternative = c("two.sided"))$p.value
  p_cont <- t.test(control0, control16, alternative = "two.sided")$p.value
  p_dec <- t.test(decline0, decline16, alternative = "two.sided")$p.value
  aa=rbind(aa, cbind(pcdh_site[pcdh_site$V4 == i, 1:3],i, 
                     mean(t(control0)), mean(t(decline0)), 
                     mean(t(control16)), mean(t(decline16)),
                     data.frame(p_0YFU), data.frame(p_16YFU),
                     data.frame(p_cont), data.frame(p_dec),
                     ann850k[i, 'UCSC_RefGene_Name']))
}
colnames(aa) <- c("chr", "pos", "pos2", "CpGsite", "mean_contr_0YFU", "mean_dec_0YFU",
                  "mean_contr_16YFU", "mean_dec_16YFU",
                  "p_0YFU", "p_16YFU", "p_control", "p_decline", "Gene")
head(aa, n=10)

write.table(aa, file = "cognitive.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)


