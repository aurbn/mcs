rm(list=ls())
source("common.R")

data.file1 = "./data/moss1_data.txt"
#samples.file1 = "./data/moss1_samples.txt"
data.file2 = "./data/moss2_data.txt"
#samples.file2 = "./data/moss2_samples.txt"
samples.fileA = "./data/mossA_samples.txt"
heatmap.file = "./results/mossA_heatmap.png"
result.file = "./results/mossA_results.txt"
out.data.file="./results/DATA_ALL_A.txt"
names.file = "./names.txt"
PV_REQ = 0.05 #Requested p-value for metabolite selection

#First set of data
data1 = read.csv(data.file1, header=T, stringsAsFactors=F)
data1[is.na(data1)] <- 0 # empty cell => no metabolite
colnames(data1) = gsub("X", "X1_", colnames(data1))
#samples1 = read.csv(samples.file1, header=T, stringsAsFactors=F)
#samples1$sample = paste0("1_", samples1$sample)
#samples1[samples1$state=="K", "state"] <- "KX"
#samples1$ctrl <- "KX"


#second set of data
data2 = read.csv(data.file2, header=T, stringsAsFactors=F)
data2[is.na(data2)] <- 0 # empty cell => no metabolite
colnames(data2) = gsub("X", "X2_", colnames(data2))
#samples2 = read.csv(samples.file2, header=T, stringsAsFactors=F)
#samples2$sample = paste0("2_", samples2$sample)
#samples2[samples2$state == "K", "state"] <- "KY"
#samples2$ctrl <- "KY"

samples = read.csv(samples.fileA, header=T, stringsAsFactors=F)
#samples2$sample = paste0("2_", samples2$sample)

#samples = rbind(samples1, samples2)
data = merge(x=data1, by.x="Name", y=data2, by.y="Name", all=TRUE)
data[is.na(data)] <- 0 #if no metabolite => C = 0

dm.m <- process_metabolome(data, samples, 10, names.file, result.file, out.data.file)

draw_heatmap(dm.m, expression(log[10](fc)), heatmap.file)

