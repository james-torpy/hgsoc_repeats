# load packages needed:
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
sampleType <- ""
types <- c("c", "t")

# define directories:
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/R/", methodName)
plotDir <- paste0(inDir, "/plots/compBarplot/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs and QC check ###

# load cCounts into list:
Counts <- list(readRDS(file=paste0(RobjectDir, "/c_RepeatCountDFs/all_cRepeatCountsDF.rds")))
names(Counts)[1] <- "cCounts"

# load tCounts and split into groups:
tCounts <- readRDS(file=paste0(RobjectDir, "/t_RepeatCountDFs/all_tRepeatCountsDF.rds"))
groups <- unique(gsub("\\/.*|\\?", "", rownames(tCounts)))
i=2
for (g in groups) {
  Counts[[i]] <- tCounts[grep(g, rownames(tCounts)),]
  names(Counts)[i] <- paste0(g, "Counts")
  i=i+1
}

# select elements of counts elements with 1 or less repeat categories:
j=1
n=1
for (i in 1:length(Counts)) {
  if (nrow(Counts[[i]]) < 3) {
    if (j==1) {
      tempCounts <- list(Counts[[i]])
      j=j+1
    } else {
      tempCounts[[j]] <- Counts[[i]]
      j=j+1
    }
  } else {
    if (n==1) {
      newCounts <- list(Counts[[i]])
      names(newCounts)[n] <- names(Counts)[i]
      n=n+1
    } else {
      newCounts[[n]] <- Counts[[i]]
      names(newCounts)[n] <- names(Counts)[i]
      n=n+1
    }
  }
}

# convert nCounts to df and add as Counts list element:
newCounts[[length(newCounts)+1]] <- do.call("rbind", newCounts)
names(newCounts)[length(newCounts)] <- "other"


### 2. Calculate percentages and plot both dfs ###

for (i in 1:length(newCounts)) {
  # calculate percentages as new data frame:
  perCounts <- apply(newCounts[[i]], 2, function(x) {
    return((x/sum(x))*100)
  })
    
  # convert df to long format for plotting:
  plotCounts <- melt(perCounts, varnames = c("type", "sample"), value.name = "percentage")
    
  # define plot colour scheme:
  #tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
  #pal(tol21rainbow)
    
  # plot data as barplot:
  if (file.exists(paste0(plotDir, "/", names(newCounts)[i], "_", sampleType, "_repeat_percent_compBarplot", ".pdf"))) {
    print(paste0(plotDir, "/", names(newCounts)[i], "_", sampleType, "_repeat_percent_compBarplot", ".pdf already exists"))
  } else {
    p <- ggplot(plotCounts, aes(x=sample, y=percentage))
    p <- p + geom_bar(stat="identity", aes(fill=type))
    p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
    pdf(file = paste0(plotDir, "/", names(newCounts)[i], "_", sampleType, "_repeat_percent_compBarplot", ".pdf"), height = 10, width = 10)
    print(p)
    dev.off()
  }
}
  

### 3. make custom repeat comp plots ###

altCounts <- list(newCounts[[7]][!rownames(newCounts[[7]]) %in% c("SINE/Alu", "SINE/MIR"),])
names(altCounts)[1] <- "altSINE"

altCounts[[2]] <- newCounts[[5]][!rownames(newCounts[[5]]) %in% "rRNA",]
names(altCounts)[2] <- "altRNA"

altCounts[[3]] <- newCounts[[2]][!rownames(newCounts[[2]]) %in% c("DNA/hAT-Charlie", "DNA/TcMar?"), ]

# halve altDNA df:
altCounts[[4]] <- altCounts[[3]][1:round(nrow(altCounts[[3]])/2),]
names(altCounts)[4] <- "altDNA2"
altCounts[[3]] <- altCounts[[3]][round(nrow(altCounts[[3]])/2):nrow(altCounts[[3]]), ]
names(altCounts)[3] <- "altDNA1"

altCounts[[5]] <- newCounts[[8]][!rownames(newCounts[[8]]) %in% c("dust", "trf"), ]
names(altCounts)[5] <- "altOther"

altCounts[[6]] <- newCounts[[4]][!rownames(newCounts[[4]]) %in% c("LTR/ERV1", "LTR/ERVL", "LTR/ERV-MaLR"), ]
names(altCounts)[6] <- "altLTR"

altCounts[[7]] <- newCounts[[3]][!rownames(newCounts[[3]]) %in% c("LINE/L1", "LINE/L2"), ]
names(altCounts)[7] <- "altLINE"

for (i in 1:length(altCounts)) {
  # calculate percentages as new data frame:
  perCounts <- apply(altCounts[[i]], 2, function(x) {
    return((x/sum(x))*100)
  })
  
  # convert df to long format for plotting:
  plotCounts <- melt(perCounts, varnames = c("type", "sample"), value.name = "percentage")
  
  # plot data as barplot:
  
  if (file.exists(paste0(plotDir, "/", names(newCounts)[i], "_", sampleType, "_repeat_percent_compBarplot", ".pdf"))) {
    print(paste0(plotDir, "/", names(altCounts)[i], "_", sampleType, "_repeat_percent_compBarplot", ".pdf already exists"))
  } else {
  p <- ggplot(plotCounts, aes(x=sample, y=percentage))
  p <- p + geom_bar(stat="identity", aes(fill=type))
  p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
  pdf(file = paste0(plotDir, "/", names(altCounts)[i], "_", sampleType, "_repeat_percent_compBarplot", ".pdf"), height = 10, width = 10)
  print(p)
  dev.off()
  }
}


### 4. make custom repeat comp plots round 2 ###

candidates <- as.character(read.csv(file=paste0(plotDir, "/repeatComp_candidates.csv"))[,1])

allRepeats <- do.call("rbind", newCounts)
allRepeats <- allRepeats[!grepl("cCounts", rownames(allRepeats)), ]
allRepeats_types <- gsub("^.*\\.", "", rownames(allRepeats))

custom2 <- allRepeats[allRepeats_types %in% candidates, ]
custom2DF <- as.data.frame(sapply(custom2, as.numeric))
rownames(custom2DF) <- rownames(custom2)

# fetch missing candidates:
missing <- gsub("_", "-", as.character(read.csv(file=paste0(plotDir, "/missing_candidates.csv"))[,1]))

for (c in missing) {
  print(grep(c, rownames(tCounts)))
}

ind <- c(50, 24, 39, 30, 7, 25, 26, 27, 28, 29, 30, 56, 60, 61, 1, 26, 10)
custom2DF <- rbind(custom2DF, tCounts[ind,])


# load in library sizes:
libFiles <- grep("subset", list.files(RobjectDir, pattern="libSize", full.names = T), value = T, invert = T)
for (i in 1:length(libFiles)) {
  if (i==1) {
    libSizes <- c(readRDS(file=libFiles[i]))
  } else {
    libSizes[i] <- readRDS(file=libFiles[i])
  }
print(i)
i=i+1
}
  
custom2CPM <- as.data.frame(t(t(custom2DF)/libSizes)*1000000)
custom2CPM$type <- gsub("^.*\\.", "", rownames(custom2CPM))
custom2CPM <- aggregate(.~type, custom2CPM, mean)
rownames(custom2CPM) <- custom2CPM$type
custom2CPM <- subset(custom2CPM, select=-type)

# fetch missing candidates:
missing <- gsub("_", "-", as.character(read.csv(file=paste0(plotDir, "/missing_candidates.csv"))[,1]))

for (c in missing) {
  print(grep(c, rownames(tCounts)))
}

ind <- c(50, 24, 39, 30, 7, 25, 26, 27, 28, 29, 30, 56, 60, 61, 1, 26, 10)
custom2Counts3 <- tCounts[ind,]
custom2CPM3 <- as.data.frame(t(t(custom2Counts3)/libSizes)*1000000)

custom2CPM <- rbind(custom2CPM, custom2CPM3)
# order by row means:
custom2CPM <- custom2CPM[order(rowMeans(custom2CPM)), ]



for (i in seq(1, 35, 7)) {
  df <- custom2CPM[i:(i+6),]
  df$Type <- factor(rownames(df))
  pDF <- melt(df, variable.name = "sample", value.name = "CPM")
  p <- ggplot(pDF, aes(x=sample, y=CPM, group=Type, colour=Type))
  p <- p + geom_line()
  p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
  print(p)
}




#j=1
#for (i in seq(1, (nrow(custom2CPM)-1), 3)) {
#  df <- as.data.frame(custom2CPM[i:(i+2), ])
#  df$type <- rownames(df)
#  pDF <- melt(df, variable.name = "sample", value.name = "CPM")
#  p <- ggplot(pDF, aes(x=sample, y=CPM))
#  p <- p + geom_bar(stat="identity", aes(fill=type), width=1, position = "dodge")
#  p <- p + theme(axis.text.x = element_text(angle = 90))
#  pdf(file=paste0(plotDir, "/custom2_CPM", j, ".pdf"))
#  print(p)
#  dev.off()
#  j=j+1
#}

###### mistake - made percentage of candidates of all repeat reads, but these values are too low.
# Calculate CPMs of these instead (with library sizes) for more appropriate values to compare



#totalCounts <- apply(allRepeats, 2, sum)
#storage.mode(totalCounts) <- "numeric"

#custom2Per <- t(t(custom2DF)/totalCounts)*100
#rownames(custom2Per) <- gsub("^.*\\.", "", rownames(custom2Per))

#for (i in seq(1, (nrow(custom2Per)-1), 2)) {
#  df <- as.data.frame(custom2Per[i:(i+1), ])
#  df$type <- rownames(df)
#  pDF <- melt(df, variable.name = "sample", value.name = "percentage")
#  p <- ggplot(pDF, aes(x=sample, y=percentage))
#  p <- p + geom_bar(stat="identity", aes(fill=type), width=0.5, position = "dodge")
#}

#cInd <- list()
#i=1
#for (c in candidates) {
#  print(c)
#  for (n in newCounts) {
#    ind <- grep(c, rownames(n))
#    if (!length(ind)==0) {
#      cInd[[i]] <- ind
#      names(ind)[i] <- names(n)[i]
#    }
#  }
#  i=i+1
#  print(i)
#}




######
# remove null values from list of count elements:
#oneCounts <- Filter(Negate(is.null), oneCounts)

# select elements of counts elements with 2 or less repeat categories:
#twoCounts <- list()
#i=1
#twoCounts <- lapply(Counts, function(x) {
#  twoCounts[[i]] <- Filter(function(y) length(y)==2, x)
#  return(twoCounts)
#  i <<- i+1
#})