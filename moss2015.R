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


dir.create("results", showWarnings = FALSE)
dir.create("results/mts", showWarnings = FALSE)
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

#Normalize by mass
data <- merge(data, samples[,c("sample","mass")], by = "sample")
data$Conc. <- data$Conc./data$mass
data$mass <- NULL

data.w <- dcast(data, Name ~ sample )
data.w <- data.w[,c("Name", as.character(samples$sample))]
rownames(data.w) <- data.w$Name
data.w$Name <- NULL
x <- apply(!is.na(data.w), 1, sum)
mdrop <- names(x[x==0])
data.w <- data.w[x!=0,]
data <- data[!(data$Name %in% mdrop),]
write.table(data.w, "results/raw.data.csv", sep = '\t', row.names = TRUE)
write.table(x[order(x, decreasing = TRUE)], "results/mrtabs_ch.csv", sep = '\t', row.names = TRUE)


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

colors <- c("K00"="red", "K05" = "blue", "K30" = "green")

for (m in mts)
{
    d <- data.w[rownames(data.w) == m,]
    d <- t(d)
    d <- data.frame(rownames(d),d, stringsAsFactors = FALSE)
    names(d) <- c("sample", "val")
    #d <- d[-1,]
    d$sample <- as.numeric(d$sample)
    d$val <- as.numeric(d$val)
    rownames(d) <- NULL
    d <- merge(d, samples[,c("sample", "state", "ctrl" )], by = "sample", all.x = TRUE)
    g <- expand.grid(unique(d$state), unique(d$ctrl))
    names(g) <- c("state", "ctrl")
    d <- merge(d, g, by = c("state", "ctrl"), all.y = T)
    d[d$ctrl == 'K5', "ctrl"] <- 'K05'
    d[d$ctrl == 'K0', "ctrl"] <- 'K00'
    
    
    d <- d[order(d$state, d$ctrl),]
    k00 <- mean(na.omit(d[d$state == "K0", "val"]))
    k05 <- mean(na.omit(d[d$state == "K5", "val"]))
    k30 <- mean(na.omit(d[d$state == "K30", "val"]))
    d <- d[d$state != "K0",]
    d <- d[d$state != "K5",]
    d <- d[d$state != "K30",]
    
    p <- ggplot(data = d, aes(x=state, y=val, fill=ctrl))+
        scale_fill_manual(values = colors)+
        geom_bar(stat="identity", position=position_dodge(), colour="black") + 
        theme(axis.text.x=element_text(angle=45, hjust=1))+
        geom_hline(aes(yintercept = k00, color = "K00"), size = 1.5)+
        geom_hline(aes(yintercept = k05, color = "K05"), size = 1.5)+
        geom_hline(aes(yintercept = k30, color = "K30"), size = 1.5)
    
    ggsave(paste0("results/mts/", gsub("/", "_", m), ".png"))
        
    
}


dm.m <- process_metabolome(d1, s1, 10, names.file, result.file, req = mts)
draw_heatmap(dm.m, expression(log[10](fc)), paste0(heatmap.file, "_01.png"), 
             mlimit = SCALE_HM, title = "Day 2")

dm.m <- process_metabolome(d2, s2, 10, names.file, result.file, req = mts)
draw_heatmap(dm.m, expression(log[10](fc)), paste0(heatmap.file, "_05.png"), 
             mlimit = SCALE_HM, title = "Day 5")

dm.m <- process_metabolome(d3, s3, 10, names.file, result.file, req = mts)
draw_heatmap(dm.m, expression(log[10](fc)), paste0(heatmap.file, "_30.png"),
             mlimit = SCALE_HM, title = "Day 30")

