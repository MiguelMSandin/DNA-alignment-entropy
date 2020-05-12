library(seqinr)
library(vegan)

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("Input fasta file must be supplied. \n", call.=FALSE)
} else {
    file <- args[1]
    
    if(length(args)>1){message("WARNING: Only first argument is passing!!")}
    
    positionInfo <- function(data, verbose=TRUE){
        info <- list()
        if(verbose){cat("Analyzing positions \n")}
        for (i in 1:(nchar(as.character(data$sequence[1])))){
            ss <- substr(data[,2],i,i)
            ss <- toupper(ss)
            df <- data.frame(base=unique(ss))
            tmp <- c()
            for ( j in 1:length(unique(ss))){
                tmp[j] <- length(grep(df[j,1], ss))
            }
            df$rep <- tmp
            rownames(df) <- as.character(df[,1])
            df[,1] <- NULL
            
            dff <- subset(df, rownames(df)!="-")
            
            info$shan[i] <- vegan::diversity(t(df))                                    # Shannon entropy of position
            info$shanc[i] <- ifelse(nrow(dff)==0, NA, vegan::diversity(t(dff)))        # Shannon entropy of position removing gaps ("-")
            info$rich[i] <- length(unique(ss))                                  # Position richness
            info$richc[i] <- ifelse(nrow(dff)==0, NA, length(rownames(dff)))    # Position richness removing gaps ("-")
            info$uniq[i] <- list(unique(ss))                                    # Unique bases in position
            info$repe[i] <- list(tmp)                                           # Repetitions of the unique bases in position
            
            if(verbose){
                if((i %% (round(nchar(as.character(data$sequence[1]))*0.1,0))) == 0){message("  ", round(i/nchar(as.character(data$sequence[1]))*100, 0), "%")}
            }
        }
        
        if(verbose){cat("Creating the data frame with the information \n")}
        
        df <- data.frame(posi=c(1:length(info$shan)),
                         shan=info$shan,
                         shanc=info$shanc,
                         rich=info$rich,
                         richc=info$richc,
                         uniq=unlist(lapply(info$uniq, function(x) paste(x, collapse="|"))),
                         repe=unlist(lapply(info$repe, function(x) paste(x, collapse="|"))))
        if(verbose){cat("Finished \n")}
        return(df)
        
    }
    
    datadf <-  seqinr::read.fasta(file, seqtype= "AA", as.string = T)
    datadf <- data.frame(name=paste(getAnnot(datadf)), sequence=paste0(datadf))
    datadf$name <- gsub(">", "", datadf$name)
    
    if(!all(duplicated(nchar(as.character(datadf$sequence)))[-1])){stop("Sequences in fasta file are not aligned.")}
    
    df <- positionInfo(datadf)
    write.table(df, paste0(gsub("\\..*{3,5}", "", file), "_positionInfo.tsv"), quote=FALSE, row.names=FALSE, sep="\t")
    cat(paste0("Data frame exported as '", gsub("\\..*{3,5}", "", file), "_positionInfo.tsv' \n"))
}
