rm(list=ls())
source("common.R")

data.file = "./data/moss1_data.txt"
samples.file = "./data/moss1_samples.txt"
heatmap.file = "./results/moss1_heatmap.png"
result.file = "./results/moss1_results.txt"
names.file = "./names.txt"
control = c('K') #list of control experiments ids
PV_REQ = 0.05 #Requested p-value for metabolite selection
MARK_K0_E0 <- ""

data = read.csv(data.file, header=T, stringsAsFactors=F)
data[is.na(data)] <- 0 # empty cell => no metabolite

samples = read.csv(samples.file, header=T, stringsAsFactors=F)

dm.m <- process_metabolome(data, samples, 2, names.file, result.file)

draw_heatmap(dm.m, expression(log[2](fc)), heatmap.file)
