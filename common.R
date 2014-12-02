process_metabolome <- function(data, samples, logbase = 2,
                               names_file = NULL, result_file=NULL,
                               out_data_file=NULL)
{
    #Collect states
    states = unique(c(as.character(samples$state))) # biological states
    all_states = as.character(states)
    controls = unique(samples[,"ctrl"])
    states = states[!(states %in% samples[,"ctrl"])] # remove control from list
    
    if(!is.null(names_file))
    {
        #convert names
        names = read.csv(names_file,sep=";", header=F, stringsAsFactors=F)
        colnames(names) <- c("tid", "hrname")
        names$tid <- gsub("(RI=.+),.+","\\1", names$tid)
        data$tid <- gsub("(RI=.+),.+","\\1", data$Name)
        data <- merge(x=data, y=names, by="tid", all.x=T, all.y=F)
        names.found <- !is.na(data$hrname)
        data[names.found, "Name"] = data[names.found, "hrname"]
        data$tid <- NULL
        data$hrname <- NULL
    }
    
    if(!is.null(out_data_file))
    {
        write.csv(data, out_data_file)
    }
       
    #Auxiliary function for statistical test
    tst <- function(row, k, e)
    {
        a = length(na.omit(row[c(k)])) 
        b = length(na.omit(row[c(e)]))
        t = list(p.value = NA) # If test cannot be performed
        try(
            {
            t <- wilcox.test(as.numeric(na.omit(row[c(k)])),
                             as.numeric(na.omit(row[c(e)])),
                             correct = F, exact = F)
            }, TRUE)
        
        if(is.finite(t$p.value))
        {
            return(t$p.value)
        }
        else
        {
            return(1)
        }   
    }
    
    Metabs <- list()  #Metabolites of interest
    
    #test each state 
    for (s in states)
    {
        control = unique(samples[samples$state == s, "ctrl"])
        sk = paste0("X", samples[samples$state %in% control,"sample"]) #coltrol samples
        se = paste0("X", samples[samples$state == s        ,"sample"]) #experiment 
        pv <- apply(data, 1, tst, sk, se)
        pv[is.na(pv)] <- 1
        Metabs <- unique(c(Metabs, as.character(data[pv < PV_REQ, "Name"])))
        
    }

    #consruct dataframe with fold changes for each found Metabolite 
    dh = data.frame(Name = data$Name) #data for heatmap
    for (s in all_states)
    {
        se = paste0("X", samples[samples$state == s, "sample"]) #experiment samples
        dh = cbind(dh, apply(data[,se], 1, mean))
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
        r[x!=0 & y==0] <- Inf #mx/mn ## OR +Inf
        r
    }
    
    #dh = as.data.frame(sweep(dh, 1, dh[,control], FUN=DF)) #normalize to control)
    #SWEEP within each control
    dh_ = NULL
    for (cn in controls)
    {
        dt_ = as.data.frame(sweep(dh[,unique(samples[samples$ctrl==cn,"state"])],
                              1, dh[,cn], FUN=DF))
        if (is.null(dh_))
        {
            dh_ = dt_
        }else
        {
            dh_ = cbind(dh_, dt_)
        }
    }
    dh = dh_
    rm(dh_)
    
    ############################################
    
    dh.raw = dh
    dh = dh[rownames(dh) %in% Metabs,]
    dh = log(dh, base=logbase)
    dh[,controls] = list(NULL)
    
    require(ggplot2)
    require(reshape2)
    require(RColorBrewer)
    
    dm = as.matrix(dh)
    dm.m = melt(dm)

    if(!is.null(result_file))
    {    
        write.csv(dh.raw, result_file)
    }
    
    return(dm.m)
}

draw_heatmap <- function(dm, label, filename=NULL, palette=NULL)
{
    #handle infinites
    t = dm$value
    t = t[is.finite(t)]
    mm = max(abs(t))
    dinf = is.infinite(dm.m$value)
    dm$value[dinf] <- mm*sign(dm$value[dinf])
    dm$inflab <- ""
    dm$inflab[dinf] <- "X"
   # dm$inflab[!dinf] <- round(dm$value[!dinf], digits=1)
    
    
    if (is.null(palette))
    {
        #COLOR palletes setut is here
        myPalette <- colorRampPalette(rev(brewer.pal(5, "RdBu")), space="rgb")
        #myPalette <- colorRampPalette(rev(c("red", "white","blue")))
    }
    mm = max(abs(dm$value))
    p <- ggplot(dm, aes(x=Var2, y=Var1, fill=value))
    p <- p + geom_tile(aes(fill = value))
    p <- p + geom_text(aes(label = inflab))
    p <- p + scale_fill_gradientn(colours = myPalette(20), 
                                  name=label, limits=c(-mm,mm))
    p <- p + xlab("Treatment") + ylab("Metabolite")
    p <- p + scale_x_discrete(expand = c(0, 0))
    p <- p + scale_y_discrete(expand = c(0, 0))
    p <- p + coord_equal()
    p <- p + theme_bw()
    p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))

    if (!is.null(filename))
    {
        ggsave(heatmap.file, scale=3)
        dev.off()
    }
    print(p)
}