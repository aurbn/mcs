rm(list=ls())
source("common.R")

data.file = "./data/moss2015_raw.txt"
samples.file = "./data/moss2015_samples.txt"
heatmap.file = "./results/moss2015.png"
result.file = "./results/moss2015.txt"
names.file = "./names.txt"
PV_REQ = 0.05 #Requested p-value for metabolite selection
require(reshape2)

# Read the data
data = read.csv(data.file, sep = '\t', header=T, stringsAsFactors=F)
data <- data[, c("FileName", "Name", "Conc.")]
data$sample <- sapply(strsplit(data$FileName, split = '-'), "[", 2)
data$sample <- sapply(strsplit(data$sample, split = '\\.'), "[", 1)
data$FileName <- NULL
data$Name <- sub("^[?[:space:]]*", "", data$Name)
d <- dcast(data, sample+Name ~ Conc.)



samples = read.csv(samples.fileA, header=T, stringsAsFactors=F)

data = merge(x=data1, by.x="Name", y=data2, by.y="Name", all=TRUE)
data[is.na(data)] <- 0 #if no metabolite => C = 0

dm.m <- process_metabolome(data, samples, 10, names.file, result.file, out.data.file)

draw_heatmap(dm.m, expression(log[10](fc)), heatmap.file)

