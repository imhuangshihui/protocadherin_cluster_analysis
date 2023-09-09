a_site <- read.table("pcdha.txt", sep = "\t", header = FALSE)
b_site <- read.table("pcdhb.txt", sep = "\t", header = FALSE)
g_site <- read.table("pcdhg.txt", sep = "\t", header = FALSE)
###########################PREPARE#############################################
pcdha_cont_0YFU <- value[rownames(value) %in% a_site$V4, conname_0YFU]
a_control0 <- c(pcdha_cont_0YFU[,1])
for(i in 2:23){
  a_control0 <- c(a_control0,pcdha_cont_0YFU[,i])
}

pcdha_cont_16YFU <- value[rownames(value) %in% a_site$V4, conname_16YFU]
a_control16 <- c(pcdha_cont_16YFU[,1])
for (i in 2:23) {
  a_control16 <- c(a_control16, pcdha_cont_16YFU[,i])
}

pcdha_dec_0YFU <- value[rownames(value) %in% a_site$V4, decname_0YFU]
a_dec0 <- c(pcdha_dec_0YFU[,1])
for (i in 2:24) {
  a_dec0 <- c(a_dec0, pcdha_dec_0YFU[,i])
}

pcdha_dec_16YFU <- value[rownames(value) %in% a_site$V4, decname_16YFU]
a_dec16 <- c(pcdha_dec_16YFU[,1])
for (i in 2:24) {
  a_dec16 <- c(a_dec16, pcdha_dec_16YFU[,i])
}
#################++++++++++++++++++++++++++###################################
pcdhb_cont_0YFU <- value[rownames(value) %in% b_site$V4, conname_0YFU]
b_control0 <- c(pcdhb_cont_0YFU[,1])
for(i in 2:23){
  b_control0 <- c(b_control0,pcdhb_cont_0YFU[,i])
}

pcdhb_cont_16YFU <- value[rownames(value) %in% b_site$V4, conname_16YFU]
b_control16 <- c(pcdhb_cont_16YFU[,1])
for (i in 2:23) {
  b_control16 <- c(b_control16, pcdhb_cont_16YFU[,i])
}

pcdhb_dec_0YFU <- value[rownames(value) %in% b_site$V4, decname_0YFU]
b_dec0 <- c(pcdhb_dec_0YFU[,1])
for (i in 2:24) {
  b_dec0 <- c(b_dec0, pcdhb_dec_0YFU[,i])
}

pcdhb_dec_16YFU <- value[rownames(value) %in% b_site$V4, decname_16YFU]
b_dec16 <- c(pcdhb_dec_16YFU[,1])
for (i in 2:24) {
  b_dec16 <- c(b_dec16, pcdhb_dec_16YFU[,i])
}
###########################+++++++++++++++++++++++++###########################
pcdhg_cont_0YFU <- value[rownames(value) %in% g_site$V4, conname_0YFU]
g_control0 <- c(pcdhg_cont_0YFU[,1])
for(i in 2:23){
  g_control0 <- c(g_control0,pcdhg_cont_0YFU[,i])
}

pcdhg_cont_16YFU <- value[rownames(value) %in% g_site$V4, conname_16YFU]
g_control16 <- c(pcdhg_cont_16YFU[,1])
for (i in 2:23) {
  g_control16 <- c(g_control16, pcdhg_cont_16YFU[,i])
}

pcdhg_dec_0YFU <- value[rownames(value) %in% g_site$V4, decname_0YFU]
g_dec0 <- c(pcdhg_dec_0YFU[,1])
for (i in 2:24) {
  g_dec0 <- c(g_dec0, pcdhg_dec_0YFU[,i])
}

pcdhg_dec_16YFU <- value[rownames(value) %in% g_site$V4, decname_16YFU]
g_dec16 <- c(pcdhg_dec_16YFU[,1])
for (i in 2:24) {
  g_dec16 <- c(g_dec16, pcdhg_dec_16YFU[,i])
}

################################################################################
plotdata_0YFU <- data.frame(
  cluster = rep(c("PCDHA", "PCDHB", "PCDHG" ), times = c(10293, 11280, 15463)),
  diagnosis = rep(c("control(n=23)", "decline(n=24)", "control(n=23)", "decline(n=24)",
                    "control(n=23)", "decline(n=24)"), 
                  times = c(5037, 5256, 5520, 5760, 7567, 7896)),
  val = c(a_control0, a_dec0, b_control0, b_dec0, g_control0, g_dec0)
)

t.test(val ~ diagnosis, data = plotdata_0YFU[plotdata_0YFU$cluster == "PCDHA",],
       alternative = "two.sided")
t.test(val ~ diagnosis, data = plotdata_0YFU[plotdata_0YFU$cluster == "PCDHB",],
       alternative = "two.sided")
t.test(val ~ diagnosis, data = plotdata_0YFU[plotdata_0YFU$cluster == "PCDHG",],
       alternative = "two.sided")

p_0YFU <- ggboxplot(plotdata_0YFU, x = "cluster", y = "val",
                    color = "diagnosis", palette = "jco", 
                    title = "DNA methylation levels of PCDH clusters between cognitive and non-cognitive patients in 0YFU") + 
  theme(plot.title = element_text(hjust = 0.5))
p_0YFU + stat_compare_means(aes(group = diagnosis, 
                                label = paste0(..method.., ", p = " ,..p.format.., "\n", ..p.signif..)), 
                            label.x = 1.5, label.y = 0.9,
                            method = "t.test")

################Do the same thing for 16YFU/Control/Decline Group###############
plotdata_16YFU <- data.frame(
  cluster = rep(c("PCDHA", "PCDHB", "PCDHG" ), times = c(10293, 11280, 15463)),
  diagnosis = rep(c("control(n=23)", "decline(n=24)", "control(n=23)", "decline(n=24)",
                    "control(n=23)", "decline(n=24)"), 
                  times = c(5037, 5256, 5520, 5760, 7567, 7896)),
  val = c(a_control16, a_dec16, b_control16, b_dec16, g_control16, g_dec16)
)

t.test(val ~ diagnosis, data = plotdata_16YFU[plotdata_16YFU$cluster == "PCDHA",],
       alternative = "two.sided")
t.test(val ~ diagnosis, data = plotdata_16YFU[plotdata_16YFU$cluster == "PCDHB",],
       alternative = "two.sided")
t.test(val ~ diagnosis, data = plotdata_16YFU[plotdata_16YFU$cluster == "PCDHG",],
       alternative = "two.sided")

p_16YFU <- ggboxplot(plotdata_16YFU, x = "cluster", y = "val",
                    color = "diagnosis", palette = "jco", 
                    title = "DNA methylation levels of PCDH clusters between cognitive and non-cognitive patients in 16YFU") +
  theme(plot.title = element_text(hjust = 0.5))
p_16YFU + stat_compare_means(aes(group = diagnosis, 
                                 label = paste0(..method.., ", p = " ,..p.format.., "\n", ..p.signif..)), 
                             label.x = 1.5, label.y = 0.9,
                             method = "t.test")
################################################################################
plotdata_cont <- data.frame(
  cluster = rep(c("PCDHA", "PCDHB", "PCDHG" ), times = c(10074, 11040, 15134)),
  year = rep(c("0YFU", "16YFU", "0YFU", "16YFU", "0YFU", "16YFU"), 
                  times = c(5037, 5037, 5520, 5520, 7567, 7567)),
  val = c(a_control0, a_control16, b_control0, b_control16, g_control0, g_control16)
)
t.test(val ~ year, data = plotdata_cont[plotdata_cont$cluster == "PCDHA",],
       alternative = "two.sided")
t.test(val ~ year, data = plotdata_cont[plotdata_cont$cluster == "PCDHB",],
       alternative = "two.sided")
t.test(val ~ year, data = plotdata_cont[plotdata_cont$cluster == "PCDHG",],
       alternative = "two.sided")

p_cont <- ggboxplot(plotdata_cont, x = "cluster", y = "val",
                    color = "year", palette = "jco", 
                    title = "DNA methylation levels of PCDH clusters in patients without cognitive impairment during age changes") +
  theme(plot.title = element_text(hjust = 0.5))
p_cont + stat_compare_means(aes(group = year, 
                                label = paste0(..method.., ", p = " ,..p.format.., "\n", ..p.signif..)), 
                            label.x = 1.5, label.y = 0.9,
                            method = "t.test")
################################################################################
plotdata_dec <- data.frame(
  cluster = rep(c("PCDHA", "PCDHB", "PCDHG" ), times = c(10512, 11520, 15792)),
  year = rep(c("0YFU", "16YFU", "0YFU", "16YFU", "0YFU", "16YFU"), 
                  times = c(5256, 5256, 5760, 5760, 7896, 7896)),
  val = c(a_dec0, a_dec16, b_dec0, b_dec16, g_dec0, g_dec16)
)

t.test(val ~ year, data = plotdata_dec[plotdata_dec$cluster == "PCDHA",],
       alternative = "two.sided")
t.test(val ~ year, data = plotdata_dec[plotdata_dec$cluster == "PCDHB",],
       alternative = "two.sided")
t.test(val ~ year, data = plotdata_dec[plotdata_dec$cluster == "PCDHG",],
       alternative = "two.sided")

p_dec <- ggboxplot(plotdata_dec, x = "cluster", y = "val",
                    color = "year", palette = "jco", 
                    title = "DNA methylation levels of PCDH clusters in patients with cognitive impairment during age changes") +
  theme(plot.title = element_text(hjust = 0.5))
p_dec + stat_compare_means(aes(group = year, 
                               label = paste0(..method.., ", p = " ,..p.format.., "\n", ..p.signif..)), 
                           label.x = 1.5, label.y = 0.9,
                           method = "t.test")

