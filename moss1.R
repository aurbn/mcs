setwd("~/mcs/Metabolome/")
data.file = "./data/moss1_data.txt"
samples.file = "./data/moss1_samples.txt"
heatmap.file = "./results/moss1_heatmap.png"
result.file = "./results/moss1_results.txt"
control = c('K') #list of control experiments ids
PV_REQ = 0.05 #Requested p-value for metabolite selection

data = read.csv(data.file, header=T)
data[is.na(data)] <- 0 # empty cell => no metabolite

samples = read.csv(samples.file, header=T)
states = unique(samples$state) # biological states
all_states = as.character(states)
states = states[!(states %in% control)] # remove control from list

#Auxiliary function for statistical test
tst <- function(row, k, e)
{
    t <- wilcox.test(as.numeric(row[c(k)]),
                     as.numeric(row[c(e)]),correct = F, exact = F)
    t$p.value
}

Metabs <- list()  #Metabolites of interest

#test each state 
for (s in states)
{
    sk = paste0("X", samples[samples$state %in% control,"sample"]) #coltrol samples
    se = paste0("X", samples[samples$state == s        ,"sample"]) #experiment samples
#     d  = data.frame(Name = data$Name, data[,sk], data[,se])
#     d$pv <- apply(d, 1, tst, sk, se) # test
#     d[is.na(d$pv), "pv"] <- 1
#     Metabs <- unique(c(Metabs, as.character(d[d$pv < PV_REQ, "Name"])))
    pv <- apply(data, 1, tst, sk, se)
    pv[is.na(pv)] <- 1
    Metabs <- unique(c(Metabs, as.character(data[pv < PV_REQ, "Name"])))
    
}

# d = data[data$Name %in% Metabs,]
# dh = data.frame(Name = d$Name) #data for heatmap
# for (s in all_states)
# {
#     se = paste0("X", samples[samples$state == s, "sample"]) #experiment samples
#     dh = cbind(dh, apply(d[,se], 1, mean))
# }
# names(dh) = c("Name", all_states)
# rownames(dh) <- dh$Name
# dh$Name <- NULL

#consruct dataframe with fold changes for each found Metabolite 
d = data[data$Name,]
dh = data.frame(Name = d$Name) #data for heatmap
for (s in all_states)
{
    se = paste0("X", samples[samples$state == s, "sample"]) #experiment samples
    dh = cbind(dh, apply(d[,se], 1, mean))
}
names(dh) = c("Name", all_states)
rownames(dh) <- dh$Name

dh$Name <- NULL

############################################
##Zero control workaround
#dh[dh[control]==0, control] <- ZEROCOUNT 
#dh <- dh+ZEROCOUNT

##Somethind bad is going here: you need to fix it!
mn = min(dh[dh>0])
mx = max(dh[dh>0])

DF <- function(x,y){
    r = x/y
    r[x==0 & y==0] <- 1
    r[x!=0 & y==0] <- mx/mn
    r
}

dh = as.data.frame(sweep(dh, 1, dh[,control], FUN=DF)) #normalize to control)
############################################

dh.raw = dh
dh = dh[rownames(dh) %in% Metabs,]
dh = log2(dh)
dh[,control] = NULL

require(ggplot2)
require(reshape2)
require(RColorBrewer)

dm = as.matrix(dh)
dm.m = melt(dm)

#handle infinites
t = dm.m$value
t = t[is.finite(t)]
mm = max(abs(t))
dinf = is.infinite(dm.m$value)
dm.m$value[dinf] <- mm*sign(dm.m$value[dinf])

#COLOR palletes setut is here
#myPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(5, "RdBu")), space="rgb")
#myPalette <- colorRampPalette(rev(c("red", "white","blue")))

p <- ggplot(dm.m, aes(x=Var2, y=Var1, fill=value))
p <- p + geom_tile()
p <- p + scale_fill_gradientn(colours = myPalette(100), name=expression(log[2](fc)),
                              limits=c(-mm,mm))
p <- p + xlab("Treatment") + ylab("Metabolite")
p <- p + scale_x_discrete(expand = c(0, 0))
p <- p + scale_y_discrete(expand = c(0, 0))
p <- p + coord_equal()
p <- p + theme_bw() 

ggsave(heatmap.file, plot=p)
dev.off()
print(p)

write.csv(dh.raw, result.file)

