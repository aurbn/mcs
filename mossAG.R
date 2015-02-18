rm(list=ls())
source("common.R")

data.file1 = "./data/moss1_data.txt"
data.file2 = "./data/moss2_data.txt"
samples.fileA = "./data/mossAG_samples.txt"
heatmap.file = "./results/mossAG_heatmap.png"
result.file = "./results/mossA_results.txt"
out.data.file="./results/DATA_ALL_AG.txt"
names.file = "./names.txt"
PV_REQ = 0.05 #Requested p-value for metabolite selection
MARK_K0_E0 <- ""

#First set of data
data1 = read.csv(data.file1, header=T, stringsAsFactors=F)
data1[is.na(data1)] <- 0 # empty cell => no metabolite
colnames(data1) = gsub("X", "X1_", colnames(data1))

#second set of data
data2 = read.csv(data.file2, header=T, stringsAsFactors=F)
data2[is.na(data2)] <- 0 # empty cell => no metabolite
colnames(data2) = gsub("X", "X2_", colnames(data2))

samples = read.csv(samples.fileA, header=T, stringsAsFactors=F)

data = merge(x=data1, by.x="Name", y=data2, by.y="Name", all=TRUE)
data[is.na(data)] <- 0 #if no metabolite => C = 0

dm.m <- process_metabolome(data, samples, 10, names.file, result.file, out.data.file)

draw_heatmap(dm.m, expression(log[10](fc)), heatmap.file)
