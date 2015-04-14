rm(list=ls())
source("common.R")
require(reshape2)

data.file = "./data/moss2015_raw.txt"
data.file.old = "./data/moss2_2015.data"

samples.file = "./data/moss2015_samples.txt"
heatmap.file = "./results/moss2015"
result.file = "./results/moss2015.txt"
names.file = "./names.txt"
PV_REQ = 0.05 #Requested p-value for metabolite selection
SCALE_HM <- 4.5
MARK_K0_E0 <- "+"


BAD_CHROMS = c(57, 58, 91, 92, 82)

# Read the data
data = read.csv(data.file, sep = '\t', header=T, stringsAsFactors=F)
data <- data[, c("FileName", "Name", "Conc.")]
data$sample <- sapply(strsplit(data$FileName, split = '-'), "[", 2)
data$sample <- sapply(strsplit(data$sample, split = '\\.'), "[", 1)
data$FileName <- NULL
data$Name <- sub("^[?[:space:]]*", "", data$Name)
full_names <- unique(data$Name)
full_samples <- unique(data$sample)
full_set <- expand.grid(Name = full_names, sample = full_samples)
data <- merge(data, full_set, by = c("Name", "sample"), all = TRUE)

samples <- read.csv(samples.file, header=T, stringsAsFactors=F)
samples <- samples[!(samples$sample %in% BAD_CHROMS),]
rownames(samples) <- samples$sample

data.w <- dcast(data, Name ~ sample )
data.w <- data.w[,c("Name", as.character(samples$sample))]
write.table(data.w, "results/raw.data.csv", sep = '\t', row.names = FALSE)

# da.na <- is.na(data.w[,-1])
# nas <- apply(da.na, 2, sum)
# nas <- nas/nrow(data.w)
# names(nas>0.7)
# bad_chroms <- unique(c(BAD_CHROMS, names(nas>0.7)))
# samples <- samples[!(samples$sample %in% bad_chroms),]

d1 <- data[data$sample %in% 1:36,]
d2 <- data[data$sample %in% 37:72,]
d3 <- data[data$sample %in% 81:106,]

d1 <- dcast(d1, Name ~ sample, value.var = "Conc.")
d2 <- dcast(d2, Name ~ sample, value.var = "Conc.")
d3 <- dcast(d3, Name ~ sample, value.var = "Conc.")

names(d1) <- paste0("X", names(d1)) # Compatibility with old code
names(d1)[1] <- "Name"
d1[is.na(d1)] <- 0 # empty cell => no metabolite
#d1[d1==""] <- 0
names(d2) <- paste0("X", names(d2)) # Compatibility with old code
names(d2)[1] <- "Name"
d2[is.na(d2)] <- 0 # empty cell => no metabolite
#d2[d2==""] <- 0
names(d3) <- paste0("X", names(d3)) # Compatibility with old code
names(d3)[1] <- "Name"
d3[is.na(d3)] <- 0 # empty cell => no metabolite
#d3[d3==""] <- 0

#d3 <- read.csv(data.file.old, header=T, stringsAsFactors=F)
#d3[is.na(d3)] <- 0 # empty cell => no metabolite


s1 <- samples[samples$ctrl == "K0",]
s1 <- s1[order(s1$state),]
s2 <- samples[samples$ctrl == "K5",]
s2 <- s2[order(s2$state),]
s3 <- samples[samples$ctrl == "K30",]
s3 <- s3[order(s3$state),]

mts <- c()

dm.m <- process_metabolome(d1, s1, 10, names.file, result.file)
mts <- c(mts, rownames(dm.m))

dm.m <- process_metabolome(d2, s2, 10, names.file, result.file)
mts <- c(mts, rownames(dm.m))

dm.m <- process_metabolome(d3, s3, 10, names.file, result.file)
mts <- c(mts, rownames(dm.m))


dm.m <- process_metabolome(d1, s1, 10, names.file, result.file, req = mts)
draw_heatmap(dm.m, expression(log[10](fc)), paste0(heatmap.file, "_01.png"), 
             mlimit = SCALE_HM, title = "Day 2")

dm.m <- process_metabolome(d2, s2, 10, names.file, result.file, req = mts)
draw_heatmap(dm.m, expression(log[10](fc)), paste0(heatmap.file, "_05.png"), 
             mlimit = SCALE_HM, title = "Day 5")

dm.m <- process_metabolome(d3, s3, 10, names.file, result.file, req = mts)
draw_heatmap(dm.m, expression(log[10](fc)), paste0(heatmap.file, "_30.png"),
             mlimit = SCALE_HM, title = "Day 30")

