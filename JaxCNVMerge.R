##########################################################################################
################################ function defination begin ###############################
is.installed <- function(requirePackage){ 
  is.element(requirePackage, installed.packages()[, 1])
}

dbscanMerge <- function(oneTypeData){
  if(nrow(oneTypeData) == 1){
    return(oneTypeData)
  }
  if(nrow(oneTypeData) == 2){
    return(twoRegionMerge(oneTypeData))
  }
  dataNum <- nrow(oneTypeData)
  mergedResults <- data.frame(chr = "", start = 0, end = 0, type = "", subtype = "")
  theLengths <- oneTypeData$end - oneTypeData$start
  theLength <- mean(oneTypeData$end - oneTypeData$start)
  theMeanDensity <- (oneTypeData$end[dataNum] - oneTypeData$start[1]) / sum(theLengths)
  theDist <- matrix(theMeanDensity + 1, nrow = dataNum, ncol = dataNum)
  theDistNeigh <- oneTypeData$start[2 : dataNum] - oneTypeData$end[1 : (dataNum - 1)]
  for(i in 1 : (dataNum - 1)){
    theDist[i, i] <- 0
    theDist[i, i + 1] <- (theLengths[i] + theLengths[i + 1] + theDistNeigh[i]) / (theLengths[i] + theLengths[i + 1])
    theDist[i + 1, i] <- theDist[i, i + 1]
  }
  theDist[dataNum, dataNum] <- 0
  theDist[which(theDist > 3, 2)] <- theMeanDensity + 1 
  theDistUsed <- as.dist(theDist)
  theRes <- dbscan(theDistUsed, minPts = 2, eps = theMeanDensity)
  
  oneCluData <- oneTypeData[which(theRes$cluster == 0), ]
  mergedResults <- rbind(mergedResults, oneCluData)
  theCluster <- unique(theRes$cluster[which(theRes$cluster != 0)])
  if(length(theCluster) > 0){
    for(clu in theCluster){
      oneCluData <- oneTypeData[which(theRes$cluster == clu), ]
      thestart = min(oneCluData$start)
      theend = max(oneCluData$end)
      mergedResults <- rbind(mergedResults, data.frame(chr = oneTypeData$chr[1], start = thestart, end = theend, 
                                                       type = oneTypeData$type[1], subtype = oneTypeData$subtype[1], 
                                                       stringsAsFactors = F))
    }
  }
  mergedResults <- mergedResults[-1, ]
  return(mergedResults)
}

bedRegionMerge <- function(oneTypeData){
  # get the distance of each pair
  oneTypeData <- oneTypeData[order(oneTypeData$start), ]
  if(nrow(oneTypeData) == 1){
    return(oneTypeData)
  }
  dataNum <- nrow(oneTypeData)
  mergedResults <- data.frame(chr = "", start = 0, end = 0, type = "", subtype = "")
  theDistNeigh <- oneTypeData$start[2 : dataNum] - oneTypeData$end[1 : (dataNum - 1)]

  ## seperate the data into different part by dist = DistCannotMerge,
  ## then in each part, using the dbscan to merge, threshold is meandensity
  theIndexLargerCannotMerge <- which(theDistNeigh > DistCannotMerge)
  if(length(theIndexLargerCannotMerge) == 0){
    mergedResults <- rbind(mergedResults, dbscanMerge(oneTypeData))
    mergedResults <- mergedResults[-1, ]
    return(mergedResults)
  }
  firstPartIndex <- 1 : theIndexLargerCannotMerge[1]
  firstPartData   <- oneTypeData[firstPartIndex, ]
  mergedResults <- rbind(mergedResults, dbscanMerge(firstPartData))
  for(j in 2 : length(theIndexLargerCannotMerge)){
    if(j > length(theIndexLargerCannotMerge)){
      break()
    }
    partIndex <- (theIndexLargerCannotMerge[j - 1] + 1) : theIndexLargerCannotMerge[j]
    partData   <- oneTypeData[partIndex, ]
    mergedResults <- rbind(mergedResults, dbscanMerge(partData))
  }
  lastPartIndex <- (theIndexLargerCannotMerge[length(theIndexLargerCannotMerge)] + 1) : dataNum
  lastPartData  <- oneTypeData[lastPartIndex, ]
  mergedResults <- rbind(mergedResults, dbscanMerge(lastPartData))
  mergedResults <- mergedResults[-1, ]
  return(mergedResults)
}

twoRegionMerge <- function(oneTypeData){
  mergedResults <- data.frame(chr = "", start = 0, end = 0, type = "", subtype = "")
  theDist <- abs(oneTypeData$end[1] - oneTypeData$start[2])
  if(theDist > DistCannotMerge){
    return(oneTypeData)
  }
  theLength <- oneTypeData$end - oneTypeData$start
  MeanTheLength <- mean(theLength)
  theFold1 <- theDist / MeanTheLength
  theFold2 <- theDist / min(theLength)
  theFold3 <- theDist / max(theLength)
  if((theFold1 < 1 & theFold2 < 3) | theFold3 < 0.1){ 
    start = min(oneTypeData$start[1], oneTypeData$start[2])
    end = max(oneTypeData$end[1], oneTypeData$end[2])
    mergedResults <- rbind(mergedResults, data.frame(chr = oneTypeData$chr[1], start = start, end = end,
                                                     type = oneTypeData$type[1], subtype = oneTypeData$subtype[1],
                                                     stringsAsFactors =F))
  }else{
    mergedResults <- rbind(mergedResults, oneTypeData)
  }
  mergedResults <- mergedResults[-1, ]
  return(mergedResults)
}

getTheMergeRes <- function(theDataUseToMerge){
  if(nrow(theDataUseToMerge) == 1){
    return(theDataUseToMerge)
  }else{
    preNum <- nrow(theDataUseToMerge)
    oneMergedRes <- bedRegionMerge(theDataUseToMerge)
    oneMergedRes <- oneMergedRes[order(oneMergedRes$start), ]
    while(nrow(oneMergedRes) < preNum){
      preNum <- nrow(oneMergedRes)
      oneMergedRes <- bedRegionMerge(oneMergedRes)
      oneMergedRes <- oneMergedRes[order(oneMergedRes$start), ]
    }
    return(oneMergedRes)
  }
}
################################# function defination end ################################
##########################################################################################

##########################################################################################
################# check if the required packages have been installed #####################
if(!is.installed("dbscan")){ install.packages("dbscan", repos="http://cran.us.r-project.org") }
if(!is.installed("data.table")){ install.packages("data.table", repos="http://cran.us.r-project.org") }
library(dbscan)
library(data.table)
##########################################################################################

##########################################################################################
###################################### argument parsing ##################################
DistCannotMerge <- 3000000
oneBed <- ""
theHelpMessge = 
  paste(" The required packages are \"dbscan\" and \"data.table\" \n",
        " The usage of JaxCNVMerge is like: \n", 
        " \"Rscript --vanilla JaxCNVMerge.R -md 3000000 -i filename\". The output file is filename.merge.bed\n",
        " Arguments: \n",
        "  --max_distance or -md  (option), numeric, distance threshold in merging, default s 3000000 \n",
        "  --bed or -i            (required), string, the bed file of the CNV fragments \n",
        "  --help or -h,          print the help messgae \n", sep = "")
args <- commandArgs(TRUE)
if(length(args) < 1) {
  args <- c("--help")
}
## Help section
if("--help" %in% args | "-h" %in% args) {
  cat(theHelpMessge)
  q(save="no")
}
if("--max_distance" %in% args){
  argIndex <- which(args == "--max_distance")
  DistCannotMerge <- as.numeric(args[argIndex + 1])
}else if("-md" %in% args){
  argIndex <- which(args == "-md")
  DistCannotMerge <- as.numeric(args[argIndex + 1])
}
if("--bed " %in% args){
  argIndex <- which(args == "--bed")
  oneBed <- args[argIndex + 1]
}else if("-i" %in% args){
  argIndex <- which(args == "-i")
  oneBed <- args[argIndex + 1]
}else{
  cat(theHelpMessge)
  q(save="no")
}
print(paste("the input file is :", oneBed, sep = ''))
print(paste("the output file is :", oneBed, ".merge.bed", sep = ""))
print(paste("the max_distance is: ", format(DistCannotMerge, scientific = F), sep = ''))
###############################################################################################

###############################################################################################
########################### merge CNVs in Bedfile bedfile #####################################
mergedResults <- data.frame(chr = "", start = 0, end = 0, type = "", subtype = "")
theData <- fread(oneBed, header = F, sep = "\t")
colnames(theData) <- c("chr", "start", "end", "type", "subtype")
chrs <- unique(theData$chr)
for(oneChr in chrs){  ## for CNVs in a chromosome
  print(paste("-------", oneChr, sep = ""))
  oneChrData <- theData[which(theData$chr == oneChr), ]
  oneChrData <- oneChrData[order(oneChrData$start), ]
  
  if(nrow(oneChrData) == 1){
    mergedResults <- rbind(mergedResults, oneChrData)
    next()
  }
  
  currType <- oneChrData$type[1]
  currEnd  <- oneChrData$end[1]
  theDataUseToMerge <- oneChrData[1, ]
  
  for(i in 2 : nrow(oneChrData)){
    if(oneChrData$type[i] == currType & (oneChrData$start[i] - currEnd) <= DistCannotMerge){
      theDataUseToMerge <- rbind(theDataUseToMerge, oneChrData[i, ])
    }else{
      mergedResults <- rbind(mergedResults, getTheMergeRes(theDataUseToMerge))
      currType <- oneChrData$type[i]
      theDataUseToMerge <- oneChrData[i, ]
    }
    currEnd  <- oneChrData$end[i]
  }
  if(nrow(theDataUseToMerge) > 0){
    mergedResults <- rbind(mergedResults, getTheMergeRes(theDataUseToMerge))
  }
}
mergedResults <- mergedResults[-1, ]
lengths <- mergedResults$end - mergedResults$start
mergedResults <- mergedResults[which(lengths > 46000), ]  ## filter the lengths
fwrite(mergedResults, file = paste(oneBed, ".merge.bed", sep = ''), quote = F, append = F, 
       row.names = F, col.names = F, sep = '\t')
###############################################################################################


