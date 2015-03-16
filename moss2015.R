rm(list=ls())
source("common.R")

data.file = "./data/moss2015_raw.txt"
data.file.old = "./data/moss2_2015.data"

samples.file = "./data/moss2015_samples.txt"
heatmap.file = "./results/moss2015.png"
result.file = "./results/moss2015.txt"
names.file = "./names.txt"
PV_REQ = 0.05 #Requested p-value for metabolite selection
MARK_K0_E0 <- "+"
require(reshape2)

# Read the data
data = read.csv(data.file, sep = '\t', header=T, stringsAsFactors=F)
data <- data[, c("FileName", "Name", "Conc.")]
data$sample <- sapply(strsplit(data$FileName, split = '-'), "[", 2)
data$sample <- sapply(strsplit(data$sample, split = '\\.'), "[", 1)
data$FileName <- NULL
data$Name <- sub("^[?[:space:]]*", "", data$Name)
d <- dcast(data, Name ~ sample, value.var = "Conc.")

names(d) <- paste0("X", names(d)) # Compatibility with old code
names(d)[1] <- "Name"

data = read.csv(data.file.old, header=T, stringsAsFactors=F)
data[is.na(data)] <- 0 # empty cell => no metabolite

data = merge(d, data, by = "Name")
data[is.na(data)] <- 0  # It's wrong remake it!!!!



samples = read.csv(samples.file, header=T, stringsAsFactors=F)



dm.m <- process_metabolome(data, samples, 10, names.file, result.file)

draw_heatmap(dm.m, expression(log[10](fc)), heatmap.file)
