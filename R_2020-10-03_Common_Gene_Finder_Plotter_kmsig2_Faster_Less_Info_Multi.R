rm(list=ls())

#COMMON GENE FINDER

#Install all packages necessary to call these libraries
library(survminer)
library(survival)
library(MASS)
library(enrichR)
library(combinat) 
library(gtools)
library(enrichR)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorRamps)
library(colorspace)
library(plyr)
library(lubridate)
library(grid)
library(ggtext)
library(ggthemes)
library (lubridate)
library(combinat)
library(qpcR)
library(foreach)
library(doParallel)
library(doSNOW)
library(data.table)
library(feather)

wrapper <- function(x, ...) #Function that wraps text after a certain number of characters
{
  paste(strwrap(x, ...), collapse = "<br>")
}

comb = function(n, x) {
  factorial(n) / factorial(n-x) / factorial(x)
}


#What gene is the trak1 cancer data scaled to? USER_INPUT
geneOfInterest = "LMNB1"

#exact name of the cancer in your data file USER_INPUT
cancerListNames = list(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal", #IGNORE
  "Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma", #IGNORE
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung", #IGNORE
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  "Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_&_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)

#Creates list of names for file read-in's, and a list of the trak1 dataframes once read in
#Modify as necessary if your trak1 files labeled in a different manner
trak1CancerListNames = paste("trak1_", cancerListNames, "_",geneOfInterest, ".txt", sep = "")
trak1CancerList <- vector(mode = "list", length = length(cancerListNames))
maxNumberGenes = 0
kmsig2SurvivalList <- vector(mode = "list", length = length(cancerListNames))
for (i in 1:length(cancerListNames)) {
  trak1CancerList[[i]] = na.omit(read.table(trak1CancerListNames[i], header = FALSE,sep=",",stringsAsFactors = FALSE))
  if (length(trak1CancerList[[c(i,5)]]) > maxNumberGenes) {
    maxNumberGenes = length(trak1CancerList[[c(i,5)]])# + 1
  }
}



#Designate how many combinations group you want to see. 
#Example: if your trak1CancerList has 8 cancers and you just want see the shared genes in them all, just set numOfComb = 8
#Example: If you want to see all possible groups of trak1CancerList with a list of 6 cancers, type in numOfComb = 2:6
#Example: If you want to see only groups of 3 out of a list of 4 cancers, just type in numOfComb = 3, and so on. 
# USER_INPUT

numOfComb = 1:length(trak1CancerList)
#OR
# numOfComb = 1:length(trak1CancerList) #defaults to finding all combiantions of your given list


#carryPreviousGroup MUST BE A SUBSET of numOfComb.
#If you do not wish to use this function, let carryPreviousGroup = 0
#The first value of carryPreviousGroup will search for all combinations of that value
#All subsequent group sizes will base the next combination on the previous group with the largest number of shared genes
#For example, if there are 20 cancers, and numOfComb = 1:5 and carryPreviousGroup = 3:5, then the program will find all combination of groups of 20 C 1, 20 C 2, 20 C 3, and then use the group of 20 C 3 with the greatest  number of shared genes and look at which of the remaining cancers would form the next group of 4 with the greatest number of shared genes
#Speeds up computing time greatly, example: let numOfComb = 1:5 with 20 cancers, if carryPreviousGroup = 0, it takes 32.35 minutes to run
#If carryPrevious Group = 3:5, it takes 23.76 seconds and it produced the exact same answer as the old method actually when I tested it
#If carryPreviousGroup = 2:5, it takes 7.7 seconds and it still produced the exact same answer
# USER_INPUT

# carryPreviousGroup = 0
#OR for a good combo size to begin with the default numOfComb option as well
carryPreviousGroup = 4:length(trak1CancerList) #default


#Initializing dataframes and variables
genesInCommon <- data.frame(matrix(ncol = 0, nrow = maxNumberGenes))
maxGenesInCommon <- data.frame(matrix(ncol = 1, nrow = 1))
geneExponent <- data.frame(matrix(ncol = 1, nrow = 1))
cancerIndexes <- data.frame(matrix(ncol = 0, nrow = 0))
maxCancerIndexes <- data.frame(matrix(ncol = length(trak1CancerList), nrow = 1))
maxCancerRows <-  data.frame(matrix(ncol = length(numOfComb), nrow = 1))
tempDataFrame <- data.frame(matrix(ncol = 1, nrow = 1))
combIndex <- data.frame(matrix(ncol = 1, nrow = 1))
groupComb <- data.frame(matrix(ncol = length(cancerListNames), nrow = 1))
testDataFrame <- data.frame(matrix(ncol = length(cancerListNames), nrow = 1))
tempGenesInCommonList <- vector(mode = "list", length = length(numOfComb))
cShift = 0
max = 0
maxGeneShareC = 1
beforeTime = Sys.time()
maxRowAndCol = 0
totalNumComb = 0
beforeInterval = 0
beforeIntervalTime = Sys.time()
interval = 100 #the smaller the number, the faster the time remianing will update


#Estimate total number of combinations 

firstTime = TRUE;
for (i in numOfComb) {
  if (i %in% carryPreviousGroup == TRUE && firstTime == TRUE) {
    totalNumComb = totalNumComb + comb(length(cancerListNames),i)
    groupComb[i] = comb(length(cancerListNames),i)
    firstTime = FALSE
  }
  else if (i %in% carryPreviousGroup == TRUE && firstTime == FALSE) {
    totalNumComb = totalNumComb + (length(cancerListNames) - i + 1)
    groupComb[i] = (length(cancerListNames) - i + 1)
  }
  else {
    totalNumComb = totalNumComb + comb(length(cancerListNames),i)
    groupComb[i] = comb(length(cancerListNames),i)
  }
}


#setup parallel backend to use many processors
cores=detectCores()
# cl <- makeCluster(cores[1]-4) #not to overload your computer
cl <- makeCluster(7)
registerDoParallel(cl)


# groupNames <- data.frame(matrix(nrow = 1, ncol = totalNumComb))
groupNames <- data.frame(matrix(nrow = 1, ncol = 0))
maxGroupNames <- data.frame(matrix(nrow = 1, ncol = length(numOfComb)))
firstTime = TRUE;
for(i in numOfComb) { #For each combination group you specify
# foreach(i = numOfComb,.combine='cbind') %dopar% {
  geneExponent <- data.frame(matrix(ncol = 1, nrow = 1))
  
  
  #Creating cancerIndexes, which will specify which groupings of cancer the program should look at
  differenceTime <- difftime(Sys.time(), beforeTime, units='mins')
  print(c("Generating cancer indexes for group size:",i, differenceTime, "mins"))
  if (i %in% carryPreviousGroup == TRUE && firstTime == TRUE) {#if you are at the same first value in carryPreviousGroup, find all possible combos
    combIndex = combinations(length(trak1CancerList),i)
    cancerIndexes = rbind.fill(cancerIndexes, as.data.frame(combIndex))
    firstTime = FALSE
  }
  else if (i %in% carryPreviousGroup == TRUE && firstTime == FALSE) { #if you are at the later values in carryPreviousGroup, basing next group on the previous largest group
    if (dim(cancerIndexes)[2] > 1) {
      tempIndex1 = cancerIndexes[maxRowAndCol,][ , colSums(is.na(cancerIndexes[maxRowAndCol,])) == 0]
    }
    else (
      tempIndex1 = cancerIndexes[maxRowAndCol,]
    )
    l = 1
    for (j in 1: length(trak1CancerList)) { 
      if (j %in% tempIndex1 == FALSE) {
          tempIndex2 = cbind(tempIndex1, j)
          tempIndex2 = data.frame(tempIndex2)
          for (k in 1:dim(tempIndex2)[2]) {
            cancerIndexes[l + cShift,k] = tempIndex2[1,k]
            combIndex[l,k] = tempIndex2[1,k]
            # print(c("carryPreviousGroup number Not 1", i, j)) #DIAGNOSTIC TOOL
          }
          l = l +1
      }
    }
  }
  else { #If you are not on a value of carryPreviousGroup, find all combinations of groups for that value
    combIndex = combinations(length(trak1CancerList) , i)
    cancerIndexes = rbind.fill(cancerIndexes, as.data.frame(combIndex))
  }
  
  tempDataFrame1 <- data.frame(matrix(ncol = 1, nrow = maxNumberGenes))
  # for (j in (cShift + 1):dim(cancerIndexes)[1]) { #looking at just one possible group combination now
  
  if ((i - 1) %in% carryPreviousGroup == TRUE) {
    differenceTime <- difftime(Sys.time(), beforeTime, units='mins')
    print(c("Running linear (single-thread) process", i, differenceTime, "mins"))
    tempGenesInCommon <- data.frame(matrix(nrow = maxNumberGenes, ncol = 1))
    for(j in (cShift + 1):dim(cancerIndexes)[1]) {
      # print(c("Running linear (single-thread) process", i, j-cShift, dim(cancerIndexes)[1] - cShift))
      cancerSubGroupNames = ""
      sharedGeneList = trak1CancerList[[c(cancerIndexes[j,1],5)]]
      for (k in 1:dim(cancerIndexes)[2]) { #creating the column string name and finding shared genes between cancers
        cancerSubGroupNames = paste(cancerSubGroupNames, cancerListNames[cancerIndexes[j,k]], sep=", ")
        # if (i == max(numOfComb)) {
        #   cancerSubGroupNamesAbbr = paste(cancerSubGroupNamesAbbr, substring(cancerListNames[cancerIndexes[j,k]],1,2), sep="")
        # }
        sharedGeneList = Reduce(intersect, list(sharedGeneList,trak1CancerList[[c(cancerIndexes[j,k],5)]]))
      }
      cancerSubGroupNames = substring(cancerSubGroupNames,3)
      for (k in 1:length(sharedGeneList)) {
      tempGenesInCommon[k,j-cShift] = sharedGeneList[k]
      }

    names(tempGenesInCommon)[j-cShift] = cancerSubGroupNames
    }
    
    
  } else {
    if ((dim(cancerIndexes)[1] - (cShift + 1)) > 3000) {
      remainder = (dim(cancerIndexes)[1] - (cShift + 1)) %% 3000
    }
    differenceTime <- difftime(Sys.time(), beforeTime, units='mins')
    print(c("Running parallel process for group:",i,differenceTime, "mins"))
    tempGenesInCommon <- foreach(j = (cShift + 1):dim(cancerIndexes)[1], .combine = cbind) %dopar% { #looking at just one possible group combination now
      cancerSubGroupNames = ""
      # cancerSubGroupNamesAbbr = ""
      sharedGeneList = trak1CancerList[[c(cancerIndexes[j,1],5)]]
      
      for (k in 1:dim(cancerIndexes)[2]) { #creating the column string name and finding shared genes between cancers
        cancerSubGroupNames = paste(cancerSubGroupNames, cancerListNames[cancerIndexes[j,k]], sep=", ")
        # if (i == max(numOfComb)) {
        #   cancerSubGroupNamesAbbr = paste(cancerSubGroupNamesAbbr, substring(cancerListNames[cancerIndexes[j,k]],1,2), sep="")
        # }
        sharedGeneList = Reduce(intersect, list(sharedGeneList,trak1CancerList[[c(cancerIndexes[j,k],5)]]))
      }
      # cancerSubGroupNamesAbbr = substring(cancerSubGroupNamesAbbr,1)
      cancerSubGroupNames = substring(cancerSubGroupNames,3)
      
      tempDataFrame2 = data.frame((do.call(cbind.data.frame, list(sharedGeneList))))
      for (l in 1:dim(tempDataFrame2)[1]) {
        tempDataFrame1[l,1] = tempDataFrame2[l,1]
        names(tempDataFrame1)[1] = cancerSubGroupNames
      }
      return(tempDataFrame1)
  
    }
    
  }
  tempGenesInCommonList[[maxGeneShareC]] = as.data.frame(tempGenesInCommon)
  maxRowAndCol = 0
  max = 0
  
  differenceTime <- difftime(Sys.time(), beforeTime, units='mins')
  print(c("Searching for maxRowAndCol", i, differenceTime, "mins"))
  for(j in 1:dim(tempGenesInCommon)[2]) {
   if(length(na.omit(tempGenesInCommon[,j])) > max) {
     max = length(na.omit(tempGenesInCommon[,j]))
     maxRowAndCol = j + cShift
   }
    if (j %% 50000 == 0) {
      differenceTime <- difftime(Sys.time(), beforeTime, units='mins')
      print(c("Continuing to search for maxRowAndCol", j, "of", groupComb[[i]], differenceTime, "mins"))
    }
  }
  maxCancerRows[1,maxGeneShareC] = maxRowAndCol

  for (j in 1:length(na.omit(cancerIndexes[maxRowAndCol,]))) {
    maxCancerIndexes[maxGeneShareC,j] = cancerIndexes[maxRowAndCol,j]
  }

  differenceTime <- difftime(Sys.time(), beforeTime, units='mins')
  print(c("Calculating mean and SE for max groupings:",i, differenceTime, "mins"))
  for (j in 1:length(na.omit(tempGenesInCommon[,maxRowAndCol - cShift]))) { #storing the V1 exponent of trak1 for a given shared gene for each cancer in that group
    for (l in 1:length(na.omit(t(cancerIndexes[maxRowAndCol,])))) {
      rowID = which(trak1CancerList[[c(cancerIndexes[maxRowAndCol,l],5)]] == tempGenesInCommon[j,maxRowAndCol-cShift])
      geneExponent[j,l] = trak1CancerList[[c(cancerIndexes[maxRowAndCol,l],1,rowID)]]
    }
  }
  geneExponentMatrix = data.matrix(geneExponent, rownames.force = NA) #forcing the geneExponent information into a matrix form

  for(j in 1:length(na.omit(tempGenesInCommon[,maxRowAndCol-cShift]))){ #Reading the largest group into the maxGenesInCommon variable
    maxGenesInCommon[j, 3 * maxGeneShareC - 2] = tempGenesInCommon[j,maxRowAndCol-cShift]
    maxGenesInCommon[j, (3 * maxGeneShareC) - 1] = mean(geneExponentMatrix[j,], na.rm = TRUE)
    maxGenesInCommon[j, (3 * maxGeneShareC)] = sd(geneExponentMatrix[j,], na.rm = TRUE) / sqrt(length(geneExponentMatrix[j,]))
  }
  
  cShift = cShift + dim(combIndex)[1]
  max = 0
  maxGeneShareC = maxGeneShareC + 1
  combIndex <- data.frame(matrix(ncol = 1, nrow = 1))
  
}


for (i in 1:length(tempGenesInCommonList)) {
  genesInCommon <- cbind(genesInCommon, as.data.frame(tempGenesInCommonList[[i]]))
  print(c("Combining gene lists", i))
}

genesInCommon = as.data.frame(genesInCommon)


for (i in 1:length(maxGroupNames)) {
  names(maxGenesInCommon)[3* i - 2] = names(genesInCommon)[maxCancerRows[1,i]]
  names(maxGenesInCommon)[3* i - 1] = "Mean"
  names(maxGenesInCommon)[3* i - 0] = "SE"
}


afterTime = Sys.time()
differenceTime = afterTime - beforeTime
differenceTime #reports how long this section took
differenceTimeSeconds <- difftime(afterTime, beforeTime, units='secs')
processingRate = dim(cancerIndexes)[1] / as.numeric(differenceTimeSeconds)
print(c("Processing rate is", processingRate, "combinations per second"))
stopCluster(cl)

#For 32 cancers,numOfComb = 1:32, previousCarryGroup = 5:32, processing rate is 291.9 combinations/min
#It took a total of 13.88 hours, for 2,432,020 combinations

#For the latest code, same conditions, processing rate is 759.24 combinations/min. It took 5.3423 hours

#You have to have .xlsx at the end of your file name. If you refer to the same file name, you can create multiple sheets
#write.xlsx(genesInCommon, file = "Work.xlsx", sheetName = "All Genes")
#write.xlsx(maxGenesInCommon, file = "Work.xlsx", sheetName = "Max Genes", append = TRUE) USER_INPUT

# Creating txt files for both the genesInCommona and maxGenesInCommon
genesInCommonFileName = paste(toString(today()), "_genesInCommon_numOfComb_", toString(min(numOfComb)), "-", toString(max(numOfComb)), "_carryPreviousGroup_", toString(min(carryPreviousGroup)), "-", toString(max(carryPreviousGroup)),"_", geneOfInterest,sep = "")
maxGenesInCommonFileName = paste(toString(today()),"_genesInCommonMAX_numOfComb_",toString(min(numOfComb)), "-", toString(max(numOfComb)), "_carryPreviousGroup_", toString(min(carryPreviousGroup)), "-", toString(max(carryPreviousGroup)),"_", geneOfInterest,sep = "")
# 
# write.table(genesInCommon, file = genesInCommonFileName, append = FALSE, quote = TRUE, sep = " ",
#             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#             col.names = TRUE, qmethod = c("escape", "double"),
#             fileEncoding = "")

#fwrite is MUCH faster than write.table
fwrite(genesInCommon, file = genesInCommonFileName, append = FALSE, quote = TRUE, sep = " ",
       eol = "\n", na = "NA", dec = ".", row.names = TRUE,
       col.names = TRUE)

# write.table(maxGenesInCommon, file = maxGenesInCommonFileName, append = FALSE, quote = TRUE, sep = " ",
#             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#             col.names = TRUE, qmethod = c("escape", "double"),
#             fileEncoding = "")

fwrite(maxGenesInCommon, file = maxGenesInCommonFileName, append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE)

#Copy and paste those txt files into an excel and you can select groups of genes and put them onto a clipboard
#Need to write in code to allow that to work
#Alternatively just proceed and it will automatically use 
differenceTime

#-----------------------------------------------------------------------------------------------------------------
#NUMBER OF GENES V # CANCER PLOTTER

#USER_INPUT Must be a subset of numOfComb, assume at least minimum is 2
numberCancersRange = 2:max(numOfComb) 


#New dataframe initialization
genesCancerPlot <- data.frame(matrix(ncol = 3, nrow = length(trak1CancerList)))
names(genesCancerPlot)[1] = "Number of Cancers"
names(genesCancerPlot)[2] = "Added Cancer"
names(genesCancerPlot)[3] = "Number of Genes In Common"

for (i in numberCancersRange) {

  tempIndex = maxCancerIndexes[i,]
  tempIndex = t(na.omit(t(tempIndex)))
  genesCancerPlot[i,1] = length(tempIndex)
  tempName = "-";
  
  if (i == min(numberCancersRange)) {#if constructin first row of dataframe, find all cancers in that subgroup
    for (j in 1:length(tempIndex)) {
      tempName = paste(tempName, gsub("_", " ", cancerListNames[[tempIndex[j]]]), sep = ", ")
    }
    genesCancerPlot[i,2] = substring(tempName, 3, nchar(tempName))
  }
  else {
    prevTempIndex = maxCancerIndexes[i - 1,]
    prevTempIndex = t(na.omit(t(prevTempIndex)))
    j = setdiff(tempIndex, prevTempIndex)
    genesCancerPlot[i,2] = paste("+", gsub("_", " ", cancerListNames[[j]]), sep = "")
  }
  
  genesCancerPlot[i,3] = length(na.omit(maxGenesInCommon[,(i*3) - 2]))
}
genesCancerPlot = na.omit(genesCancerPlot)


genesCancerPlotName = paste(toString(today()), "_Number_Of_Genes_In_", toString(min(numberCancersRange)), "-", toString(max(numberCancersRange)), "_Cancers_", geneOfInterest, ".pdf", sep = "")

# pdf(file = genesCancerPlotName, width = 11, height = 7) #USER_INPUT. Automatically creates a PDF of the plot, comment out if you don't want it

ggplot(genesCancerPlot,aes(x = genesCancerPlot$`Number of Cancers` ,y = genesCancerPlot$`Number of Genes In Common`)) +
  geom_line(size = 0.25)+
  geom_point(shape = 18, size  = 3) +
  geom_text(aes(label=genesCancerPlot[,2]),angle = 30, hjust= -.05, vjust= -1.5, size= 2) + #If you wish to label the cancer additions, otherwise comment out
  theme_classic() +
  ggtitle(paste("Maximum Number of ", geneOfInterest, "-Scaling", " Genes Shared In Each Group of Cancer", sep = "")) +
  xlab("Number of Cancers per Group") +
  ylab("Maximum Number of Genes Shared In Cancer Group")+
  ylim(0, 1.1 * (max(genesCancerPlot[,3])))
  # theme(axis.text.x = element_text(face = "bold", color="black", size=6.5)) +
  # theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +

# dev.off() #USER_INPUT. also need to comment this out or in depending on the pdf code line as well





#----------------------------------------------------------------------------------------------------------------------------------------
#COMMON GENE PLOTTER, YOU MUST RUN THE PREVIOUS SECTION FOR THIS TO WORK

#USER_INPUT
#Write in column number of where that group is listed manually, remove"#" underenath and comment out the other one
#Must be a (multiple of 3) - 2
# columnOfInterest = 49
#OR 
# columnOfInterest = dim(maxGenesInCommon)[2] - 2 #automatically assume biggest group of in maxGenesInCommon.
#OR
limitSize = 10 #what column contains the number of shared genes you wish to see, specify in this variable

#USER_INPUT
#Write in genesInCommon or maxGenesInCommon for which dataframe you wish to look from
genesInCommonDataSet = maxGenesInCommon

#USER_INPUT, what method did you produce the kmsig2 dataframe, survfit or Cox proportional hazards?
method = "survfit"
#OR
# method = "coxph"






for (i in 1:(dim(maxGenesInCommon)[2]/3)) {
  j = (i * 3) - 2
  if (length(na.omit(maxGenesInCommon[,j])) < limitSize) {
    columnOfInterest = j - 3
    break;
  }
  else {
    columnOfInterest = j
  }
}
  
kmsig2SurvivalListNames = paste("kmsig2",cancerListNames, method, sep = "_") 
kmsig2SurvivalList <- vector(mode = "list", length = length(cancerListNames))
for (i in 1:length(cancerListNames)) {
  kmsig2SurvivalList[[i]] = read.table(kmsig2SurvivalListNames[i], header = TRUE,sep=" ",stringsAsFactors = FALSE)
}


#Generating the new dataframe genesInCommonBar where the info will be plotted
genesInCommonBar <- data.frame(matrix(ncol = 1, nrow = length(genesInCommonDataSet[,columnOfInterest])))
genesInCommonBar[,1] = genesInCommonDataSet[,columnOfInterest]
genesInCommonBar[,2] = genesInCommonDataSet[,columnOfInterest+1]
genesInCommonBar[,3] = genesInCommonDataSet[,columnOfInterest+2]
genesInCommonBar[is.na(genesInCommonBar)] <- 0 #NA replaced with 0

for(i in 1:dim(genesInCommonBar)[1]) { #any row with all 0s will be removed 
  if (genesInCommonBar[i,1] == 0 && genesInCommonBar[i,2] == 0 && genesInCommonBar[i,3] == 0) {
    genesInCommonBar <- genesInCommonBar[-c(i:dim(genesInCommonBar)[1]),]
    break
  }
}

#Generating the mean - SE and mean + SE variables
genesInCommonBar[,2] = as.numeric(genesInCommonBar[,2])
genesInCommonBar[,3] = as.numeric(genesInCommonBar[,3])
genesInCommonBar = na.omit(genesInCommonBar)
genesInCommonBar[,4] = genesInCommonBar[,2] - genesInCommonBar[,3]
genesInCommonBar[,5] = genesInCommonBar[,2] + genesInCommonBar[,3]


#arrange based on mean of exponenets, factor makes sure ggplot won't order it according to its default settings
genesInCommonBar <- genesInCommonBar %>% arrange(genesInCommonBar[,2])
genesInCommonBar[,1] <- factor(genesInCommonBar[,1], levels = genesInCommonBar[,1])

#searching for the indexes of the cancers in your selected group
for(i in 1:dim(genesInCommon)[2]) {
  if (names(genesInCommon)[i] == names(genesInCommonDataSet)[columnOfInterest]) {
    genesInCommonIndex = (i + 2) / 3
    break
  }
}

names(genesInCommonBar)[1] = names(genesInCommon)[(genesInCommonIndex * 3) - 2]
names(genesInCommonBar)[2] = "Mean_Exponent"
names(genesInCommonBar)[3] = "SE"
names(genesInCommonBar)[4] = "Lower"
names(genesInCommonBar)[5] = "Upper"

#Counting number of cancers where a gene in your selected list correlates significantly to positive or negative survival

tempIndex = cancerIndexes[genesInCommonIndex,]
tempIndex = tempIndex[ , colSums(is.na(tempIndex)) == 0]
sigSurvivalTracker <- data.frame(matrix(ncol = 2 * (dim(genesInCommonBar)[1]) + 5, nrow = length(tempIndex)))
sigSurvivalTracker[is.na(sigSurvivalTracker)] = 0
names(sigSurvivalTracker)[1] = "Cancers"
for(i in 1:dim(genesInCommonBar)[1]) {
  genesInCommonBar[i,6] = 0
  genesInCommonBar[i,7] = 0
  genesInCommonBar[i,8] = 0
  genesInCommonBar[i,9] = 0
  names(sigSurvivalTracker)[(2 * i) + 4] = paste(genesInCommonBar[i,1], "\nNeg_Surv", sep = '')
  names(sigSurvivalTracker)[(2 * i) + 5] = paste(genesInCommonBar[i,1], "\nPos_Surv", sep = '')
  # print(c("checking survival correlation for genes in common", i)) #DIAGNOSTIC_TOOL
  for(j in 1:dim(tempIndex)[2]) {
    tempDataFrame = na.omit(kmsig2SurvivalList[[c(tempIndex[,j])]])
    rowID = which(tempDataFrame[,4] == genesInCommonBar[i,1])
    sigSurvivalTracker[j,1] = cancerListNames[tempIndex[[j]]]
    if (length(rowID) > 0) { #it is possible to create rowID's of no length, we need to filter that out
      if (method == "coxph") {
        if (tempDataFrame[rowID,5] > 1 && tempDataFrame[rowID,1] < 0.05) { #negative survival in coxph files
          genesInCommonBar[i,6] = genesInCommonBar[i,6] + 1
          # print(c("Negative Sig Survival", toString(genesInCommonBar[i,1]), cancerListNames[[c(tempIndex[,j])]])) #DIAGNOSTIC_TOOL
          sigSurvivalTracker[j, (2 * i) + 4] = 1
          sigSurvivalTracker[j, 2] = sigSurvivalTracker[j, 2] + 1
        }
        else if (tempDataFrame[rowID,5] < 1 && tempDataFrame[rowID,1] < 0.05) { #positive survival in coxph files
          genesInCommonBar[i,7] = genesInCommonBar[i,7] + 1
          # print(c("Negative Sig Survival", toString(genesInCommonBar[i,1]), cancerListNames[[c(tempIndex[,j])]])) #DIAGNOSTIC_TOOL
          sigSurvivalTracker[j,(2* i) + 5] = 1
          sigSurvivalTracker[j, 3] = sigSurvivalTracker[j, 3] + 1
        }
        else { #count number of cancers the gene does not significantly predict survival
          genesInCommonBar[i,8] = genesInCommonBar[i,8] + 1
        }
      }
      else if (method == "survfit") {
        if (tempDataFrame[rowID,5] < 1 && tempDataFrame[rowID,1] < 0.05) { #negative survival in survfit files
          genesInCommonBar[i,6] = genesInCommonBar[i,6] + 1
          # print(c("Negative Sig Survival", toString(genesInCommonBar[i,1]), cancerListNames[[c(tempIndex[,j])]])) #DIAGNOSTIC_TOOL
          sigSurvivalTracker[j, (2 * i) + 4] = 1
          sigSurvivalTracker[j, 2] = sigSurvivalTracker[j, 2] + 1
        }
        else if (tempDataFrame[rowID,5] > 1 && tempDataFrame[rowID,1] < 0.05) { #positive survival in survfit files
          genesInCommonBar[i,7] = genesInCommonBar[i,7] + 1
          # print(c("Positive Sig Survival", toString(genesInCommonBar[i,1]), cancerListNames[[c(tempIndex[,j])]])) #DIAGNOSTIC_TOOL
          sigSurvivalTracker[j,(2* i) + 5] = 1
          sigSurvivalTracker[j, 3] = sigSurvivalTracker[j, 3] + 1
        }
        else { #count number of cancers the gene does not significantly predict survival
          genesInCommonBar[i,8] = genesInCommonBar[i,8] + 1
        }
      }
    }
    else if(length(rowID) == 0) { #count number of NA's, only will be seen with the survfit method
      genesInCommonBar[i,9] = genesInCommonBar[i,9] + 1
    }
  }
}
names(sigSurvivalTracker)[2] = paste("#_Neg_Surv_Genes\n", method, sep = "")
names(sigSurvivalTracker)[3] = paste("#_Pos_Surv_Genes\n", method, sep = "")
names(sigSurvivalTracker)[4] = "Color_Index"
names(sigSurvivalTracker)[5] = "Color_ID"

maxSig = max(sigSurvivalTracker[,2] + sigSurvivalTracker[,3])
colorGradient <- colorRamp(c("red", "blue"))

for(i in 1:dim(sigSurvivalTracker)[1]) {
  if (sigSurvivalTracker[i,2] == 0 && sigSurvivalTracker[i, 3] == 0) {
    sigSurvivalTracker[i,4] = -1 #allows recognition that this gene does not affect survival at all
    sigSurvivalTracker[i,5] = "black"
  }
  else {
    colIndex = ((sigSurvivalTracker[i,3] - sigSurvivalTracker[i,2]) / (2 * maxSig)) + 0.5
    sigSurvivalTracker[i,4] = colIndex
    sigSurvivalTracker[i,5] = rgb(colorGradient(sigSurvivalTracker[i,4]) , max = 255)
  }
}

names(genesInCommonBar)[6] = "Sig_Lower_Surv"
names(genesInCommonBar)[7] = "Sig_Higher_Surv"
names(genesInCommonBar)[8] = "Not_Sig"
names(genesInCommonBar)[9] = "NA_Surv"

#colorGradient is a function, an input of 0 is "red", an in put of 1 is "blue", it produces continuous colors
#For example an input of 0.5 is "purple"
colorGradient <- colorRamp(c("red", "blue"))

#Finding the highest number of cancers where a given gene affects survival
maxSig = max(genesInCommonBar$Sig_Lower_Surv + genesInCommonBar$Sig_Higher_Surv)

#creates the value from 0 to 1 for the colorGradient function
for(i in 1:dim(genesInCommonBar)[1]) {
  if (genesInCommonBar[i,6] == 0 && genesInCommonBar[i, 7] == 0) {
    genesInCommonBar[i,10] = -1 #allows recognition that this gene does not affect survival at all
    genesInCommonBar[i,11] = "black"
  }
  else {
    #colIndex = ((genesInCommonBar[i,7] - genesInCommonBar[i,6]) / (2 *(genesInCommonBar[i,6] + genesInCommonBar[i,7]))) + 0.5
    #colIndex = ((genesInCommonBar[i,7] - genesInCommonBar[i,6]) / (2 * dim(tempIndex)[2])) + 0.5
    colIndex = ((genesInCommonBar[i,7] - genesInCommonBar[i,6]) / (2 * maxSig)) + 0.5
    genesInCommonBar[i,10] = colIndex
    genesInCommonBar[i,11] = rgb(colorGradient(genesInCommonBar[i,10]) , max = 255)
  }
}


#generates point labels based on number of negatively effects groups and positivley effected groups
#Comment out points you don't want to see

if (method == "coxph") {
  for(i in 1:dim(genesInCommonBar)[1]) {
    genesInCommonBar[i,12] = paste(
      genesInCommonBar[i,6], #label # Predict Neg Survival
      genesInCommonBar[i,7], #label # Predict Pos Survival
      genesInCommonBar[i,8], #label # Not Sig Predictor
      # genesInCommonBar[i,9], #label # NA's
      sep = ",", collapse = NULL)
  }
} else if (method == "survfit") {
  for(i in 1:dim(genesInCommonBar)[1]) {
    genesInCommonBar[i,12] = paste(
    genesInCommonBar[i,6], #label # Predict Neg Survival
    genesInCommonBar[i,7], #label # Predict Pos Survival
    genesInCommonBar[i,8], #label # Not Sig Predictor
    genesInCommonBar[i,9], #label # NA's
    sep = ",", collapse = NULL)
  }
}

names(genesInCommonBar)[10] = "Color_Index"
names(genesInCommonBar)[11] = "Color_ID"
names(genesInCommonBar)[12] = "negSurv#,posSurvival#, notSig#, NA#"

minColorIndex = min(genesInCommonBar[,10])
meanColorIndex = mean(range(genesInCommonBar[,10]))
maxColorIndex = max(genesInCommonBar[,10])



wrapAroundValue = 90 #how many characters until text wraps around USER_INPUT
fontSize = "9" #Font size of title, adjust this and wrapAroundValue until it all fits USER_INPUT


genesInCommonBarTitle = paste("**Scaling gene** = ", geneOfInterest,";  ", "**Method** = ", method,
                              ";  ", "<br>Point labels = (# Predict Neg Surv, # Predict Pos Surv, # Not Sig Predictors, # NA's [if present])",
                              "<br><br>**Cancers**: ", sep = "")

#Generates markdown to read in title with colored words and different formatting
characterCount = 0
for (i in 1:dim(sigSurvivalTracker)[1]) {
  genesInCommonBarTitle = paste(genesInCommonBarTitle, "<span style='color:", toString(sigSurvivalTracker[i,5]), ";'>", 
                                gsub("_", " ", toString(sigSurvivalTracker[i,1])), sep = "") #does this work <br> </span>
  characterCount = characterCount + nchar(toString(sigSurvivalTracker[i,1]));
  if (sigSurvivalTracker[i,4] != -1) {
    genesInCommonBarTitle = paste(genesInCommonBarTitle, " (", toString(sigSurvivalTracker[i,2]), ",", toString(sigSurvivalTracker[i,3]), ")", sep = "" )
  }
  genesInCommonBarTitle = paste(genesInCommonBarTitle, "</span>, ", sep = "")
  
  if (characterCount > wrapAroundValue) {
    genesInCommonBarTitle = paste(genesInCommonBarTitle, "<br>", sep = "")
    characterCount = 0
  }
}
genesInCommonBarTitle = substring(genesInCommonBarTitle,1, nchar(genesInCommonBarTitle)-2)
genesInCommonBarTitle = paste("<span style='font-size:", fontSize, "pt'>", genesInCommonBarTitle, "</span>", sep = "" )


genesInCommonGraphName = paste(toString(today()), "_Shared_Genes_", toString(length(tempIndex)), "_Cancers_", geneOfInterest, "_", method, ".pdf", sep = "")

#USER_INPUT, read this if you want to create a PDF automatically
# pdf(file = genesInCommonGraphName, width = 11, height = 7)


if (dim(genesInCommonBar)[1] > 1) {
  #Plots the graph
  ggplot(genesInCommonBar,aes(x = genesInCommonBar[,1] ,y = Mean_Exponent, fill = genesInCommonBar[,10])) +
    scale_fill_gradient2(name = "% Predict Neg Surv,\n% Predict Pos Surv\n",
                         breaks = c(minColorIndex,meanColorIndex,maxColorIndex),
                         label = c("0%, 100%", "50%, 50%", "100%, 0%"),
                         midpoint = meanColorIndex, low = '#0000FF', mid = '#7F007F', high = '#FF0000') +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), size  = 0.5, width = 0.3, colour = genesInCommonBar[,11]) +
    geom_point(shape = 18, size  = 3, colour = genesInCommonBar[,11]) +
    geom_text(aes(label=genesInCommonBar[,12]),hjust= -.2, vjust= -.5, size= 2) +
    theme_bw() +
    theme(axis.text.x = element_text(face = "bold", color="black", size=6.5)) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    ylab("Exponent") +
    xlab("Genes") +
    # ggtitle(paste("Scaling gene = ", geneOfInterest,"     ", "Method = ", method,
                  # "               ", "Point labels = (# Predict Neg Surv, # Predict Pos Surv, # Not Sig Predictors, # NA's)",
                  # "\n\nCancers: ",
                  # wrapper(names(genesInCommonBar)[1], width = 150), sep = "")) +
    # theme(axis.title= element_text(face  = "bold")) +
    labs(title = genesInCommonBarTitle) +
    theme(plot.title = element_markdown()) 
} else {
  #Plots the graph
  ggplot(genesInCommonBar,aes(x = genesInCommonBar[,1] ,y = Mean_Exponent)) +
    xlab("Genes") +
    # scale_fill_gradient2(name = "% Predict Neg Surv,\n% Predict Pos Surv\n")
    geom_errorbar(aes(ymin = Lower, ymax = Upper), size  = 0.5, width = 0.3, colour = genesInCommonBar[,11]) +
    geom_point(shape = 18, size  = 3, colour = genesInCommonBar[,11]) +
    geom_text(aes(label=genesInCommonBar[,12]),hjust= -.2, vjust= -.5, size= 2) +
    theme_bw() +
    theme(axis.text.x = element_text(face = "bold", color="black", size=6.5)) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    labs(title = genesInCommonBarTitle) +
    theme(plot.title = element_markdown()) +
    ylab("Exponent")
}
  
# dev.off() #USER_INPUT uncomment this if the pdf function is uncommented as well




#-----------------------------------------------------------------------------------------------------
#kmsig2 FILE CREATOR USING EITHER SURVFIT OR COXPH METHODOLOGY

# rm(list=ls()) #remove any previous variables
#setwd("/Applications/UPenn/Summer2020/Pancreatic")

library(survminer)
library(survival)
library(MASS)
library(enrichR)
library(combinat)
library(gtools)
library(enrichR)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorRamps)
library(colorspace)
library(plyr)

#exact name of the cancer in your downloaded data file from the TCGA site USER_INPUT
cancerListNames = list(
  #"Acute_Myeloid_Leukemia" #Problematic, no primary tumor site
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  "Colon_and_Rectal",
  "Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  "Lower_Grade_Glioma_and_Glioblastoma",
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  "Lung",
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  "Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_&_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)

#USER_INPUT, what method did you produce the kmsig2 variable, survfit or Cox proportional hazards?
method = "coxph"
#OR
# method = "survfit"



dataSurvivalListNames = paste("Data_", cancerListNames, "_survival", ".txt", sep = "") 
dataHiSeqV2ListNames = paste("Data_",cancerListNames,"_HiSeqV2", sep = "") 
dataSurvivalList <- vector(mode = "list", length = length(cancerListNames))
dataHiSeqV2List <- vector(mode = "list", length = length(cancerListNames))
fileNameList = paste("kmsig2_", cancerListNames, "_", method, sep = "") 
for (i in 1:length(cancerListNames)) {
  print(c("reading in files", i))
  dataSurvivalList[[i]] = read.table(dataSurvivalListNames[i], header = TRUE,sep="\t",stringsAsFactors = FALSE)
  dataHiSeqV2List[[i]] = read.table(dataHiSeqV2ListNames[i], header = FALSE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
}

for (x in 1:length(cancerListNames)) {
  df.input = dataSurvivalList[[x]] #Phenotype data like survival years etc
  df.map = dataHiSeqV2List[[x]] #Matrix of genes vs patient id
  df.map <- na.omit(df.map) #remove NA if any
  ## narrowing down to only primary tumor sites
  Trk1=df.map[2:dim(df.map)[1],2:3] #data frame storing mRNA seq values only
  name=df.input[1] #patient ids from the survival file
  namepro=df.map[1,2:dim(df.map)[2]] #patient ids from mrna seq file
  patientid=substr(namepro[1],1,12) #pure patient id
  k=1
  
  for(i in 1:dim(namepro)[2])
  {
    if(substr(namepro[1,i],14,15)=='01')
    {
      print(c("finding primary tumor sites", i))
      Trk1[,k]=as.double(df.map[2:dim(df.map)[1],i+1]) #mrna reads
      patientid[k]=substr(namepro[i],1,12)
      k=k+1
    }
  }
  
  ##matching RNA seq patient data to survival data
  Trk <- data.frame(matrix(ncol = 3, nrow = 2))
  for(i in 1:length(patientid))
  {
    for(j in 1:dim(df.input)[1])
    {
      print(c("matching RNA seq", i,j))
      if(patientid[i]==df.input$X_PATIENT[j] )
      {
        Trk[i,3]=df.input$OS[j]
        Trk[i,2]=df.input$OS.time[j]/365
        break
      }
    }
  }
  
  chk=0.75*dim(Trk1)[2] 
  kmsig <- data.frame(matrix(ncol = 3, nrow = 2))
  pro_name=df.map[1,1]
  k=1
  
  if (method == "coxph") {
    for (i in 1:dim(Trk1)[1])
    {
      print(c("coxph survival", x, i))
      if(sum(Trk1[i,]==0)<chk)
      {
        Trk[,1]=as.double(t(Trk1[i,]))
        Trk[,4]=NULL
        Trk[,4] <- Trk[,1]> median(Trk[,1])
        coxph.fit <- coxph(Surv(Trk[,2],Trk[,3]) ~ Trk[,4], data = Trk)
        
        kmsig[k,1]=coef(summary(coxph.fit))[,5] # p value, specifically Pr(>|z|)
        kmsig[k,2]=coef(summary(coxph.fit))[,3] # se(coef)
        kmsig[k,3]=coef(summary(coxph.fit))[,4] #z
        kmsig[k,4]=df.map[i+1,1]
        kmsig[k,5]=coef(summary(coxph.fit))[,2] #hazard ratio
        pro_name[k]=df.map[i+1,1]
        k=k+1
      }
    }
    kmsig2 = kmsig
  
    names(kmsig2)[1] = "p-value_Coeff"
    names(kmsig2)[2] = "SE_Coef"
    names(kmsig2)[3] = "Z"
    names(kmsig2)[4] = "Gene"
    names(kmsig2)[5] = "Hazard_Ratio"
  }
  else if (method == "survfit") {
    for (i in 1:dim(Trk1)[1])
    {
      ##TEST
      print(c("survfit survival", x, i))
      
      if(sum(Trk1[i,]==0)<chk)
      {
        Trk[,1]=as.double(t(Trk1[i,]))
        Trk[,4]=NULL
        Trk[,4] <- Trk[,1]> median(Trk[,1])
        fit <- survfit(Surv(Trk[,2],Trk[,3]) ~ Trk[,4],
                       data = Trk)
        
        kmsig[k,1]=surv_pvalue(fit)$pval #p value
        kmsig[k,2]=surv_median(fit)$median[1] #median survival of patients with lower expression of that gene
        kmsig[k,3]=surv_median(fit)$median[2] #median survival of patients with higher expression of that gene
        pro_name[k]=df.map[i+1,1]
        k=k+1
      }
    }
    
    kmsig1=kmsig
    kmsig2=cbind(kmsig1,pro_name)
    kmsig2[,5]=kmsig2[,3]/kmsig2[,2] #ratio of survival time of higher expressors:lower expressors
    names(kmsig2)[1] = "p-value_Coeff"
    names(kmsig2)[2] = "med_surv_lower_expressed"
    names(kmsig2)[3] = "med_surv_higher_expressed"
    names(kmsig2)[4] = "Gene"
    names(kmsig2)[5] = "Surv_Ratio_(Col3/Col2)"
  }
  
  write.table(kmsig2, file = fileNameList[x], append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = TRUE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "") #Change Tumor Name here
}






#------------------------------------------------------------------------------------------------------
#FIND PERCENT OF STRONGLY SCALING GENES THAT PREDICT POSITIVE AND NEGATIVE SURVIVAL

kmsig2=read.table("kmsig2_Esophagus_Survival", header = TRUE,sep=" ",stringsAsFactors = FALSE)
df.strong_scaling = read.table("trak1_Esophagus_ACTA2.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)

LMNB1_km=subset(kmsig2,kmsig2$pro_name %in% df.strong_scaling[,5])

#how many genes that scale strongly predict significant survival
scaling_survival1=subset(LMNB1_km, LMNB1_km$X1<0.05 & LMNB1_km$V5<1)
scaling_survival2=subset(LMNB1_km, LMNB1_km$X1<0.05 & LMNB1_km$V5>1)


#scaling_survival1=subset(LMNB1_km, LMNB1_km$X1<0.05)
#what isthe percentage with respect to all the scaling genes
per=dim(scaling_survival1)[1]/dim(LMNB1_km)[1]
per1=dim(scaling_survival2)[1]/dim(LMNB1_km)[1]
plot(-log2(as.double(kmsig2[,1])),log2(as.double(kmsig2[,5])),pch = 19,cex=0.3,col="grey",xlab='-log2(p-value)',ylab='log2(median survival factor change)',title(c(signif(per*100,2)," % show poor survival",per1*100," % show pro survival")))

points(-log2(LMNB1_km[,1]),log2(LMNB1_km[,5]),pch = 19,cex=0.5,col="red")
points(c(-log2(0.05),-log2(0.05),-log2(0.05)),-1:1,pch = 3,cex=1)

plot(kmsig2[,1],kmsig2[,5],pch = 19,cex=0.3,col="grey",xlab='-log2(p-value)',ylab='log2(median survival factor change)',title(c(signif(per*100,2)," % show poor survival",per1*100," % show pro survival")))











#-------------------------------------------------------------------------------------------------------------------------
# INDIVIDUAL KM PLOTTER FOR A SINGLE CANCER AND SINGLE GENE

#USER_INPUT, read in your particular cancer files from TCGA
df.input = read.table("Data_Breast_survival.txt", header = TRUE,sep="\t",stringsAsFactors = FALSE) #Phenotype data like survival years etc
df.map = read.table("Data_Breast_HiSeqV2", header = FALSE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA")) #Matrix of genes vs patient id
df.map <- na.omit(df.map)

## narrowing down to only primary tumor sites
Trk1=df.map[2:dim(df.map)[1],2:3]
name=df.input[1]
namepro=df.map[1,2:dim(df.map)[2]]
patientid=substr(namepro[1],1,12)
k=1

for(i in 1:dim(namepro)[2])
{
  if(substr(namepro[1,i],14,15)=='01')
  {
    Trk1[,k]=as.double(df.map[2:dim(df.map)[1],i+1])
    patientid[k]=substr(namepro[i],1,12)
    k=k+1
    
    ##TEST
    print(c("narrowing to only primary tumor sites", k))
  }
}
j = 0
##matching RNA seq patient data to survival data
Trk <- data.frame(matrix(ncol = 3, nrow = 2))
for(i in 1:length(patientid))
{
  for(j in 1:dim(df.input)[1])
  {
    ##TEST
    print(c("matching RNA seq", i,j))
    
    if(patientid[i]==df.input$X_PATIENT[j] )
    {
      Trk[i,3]=df.input$OS[j]
      Trk[i,2]=df.input$OS.time[j]/365
      break
    }
  }
}

# searching for the gene of interest fr sideways volcano plot
gene = "SKA3"

for (j in 2:dim(df.map)[1])
{
  if(df.map[j,1]== gene) #USER_INPUT
  {
    iden=j-1
    break
  }
}


Trk[,1]=as.double(t(Trk1[iden,]))
Trk[,4]=NULL
Trk[,4] <- Trk[,1]> median(Trk[,1])
fit <- survfit(Surv(Trk[,2],Trk[,3]) ~ Trk[,4],
               data = Trk)
ggsurvplot(fit, data = Trk, risk.table = FALSE,pval = TRUE, title = gene, xlab = "Survival time in years",font.x=15,font.y=15,font.tickslab=15)








#---------------Plotting Gene Ontology enrichment-------------------------------------

library(ggrepel)
set.seed(42)
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018","KEGG_2019_Human","WikiPathways_2019_Human")
enriched <- enrichr(df.strong_scaling$V5, dbs)
plot(df.strong_scaling[,3],df.strong_scaling[,1],pch=19,xlab=parse(text="R^2"),ylab="Exponent")

an1=enriched[["GO_Biological_Process_2018"]]$Genes[enriched[["GO_Biological_Process_2018"]]$Term=="extracellular matrix organization (GO:0030198)"]
an2=enriched[["GO_Biological_Process_2018"]]$Genes[enriched[["GO_Biological_Process_2018"]]$Term=="regulation of angiogenesis (GO:0045765)"]
an3=enriched[["KEGG_2019_Human"]]$Genes[enriched[["KEGG_2019_Human"]]$Term=="Focal adhesion"]

an1s=unlist(strsplit(an1, ";"))
an1p=subset(df.strong_scaling,df.strong_scaling$V5 %in% an1s)
points(an1p[,3],an1p[,1],pch=19,col="red")
text(an1p[,3],an1p[,1],an1p[,5],pos=3)

an2s=unlist(strsplit(an2, ";"))
an2p=subset(df.strong_scaling,df.strong_scaling$V5 %in% an2s)
points(an2p[,3],an2p[,1],pch=19,col="blue")

an3s=unlist(strsplit(an3, ";"))
an3p=subset(df.strong_scaling,df.strong_scaling$V5 %in% an3s)
points(an3p[,3],an3p[,1],pch=19,col="green")
#----alternate repel text label plotting
ggplot(df.strong_scaling,aes(V3,V1,label=V5))+geom_point(color = "blue")+geom_label_repel()
enriched[["KEGG_2019_Human"]]







data <- data.frame(
  group = rep(c('affluence', 'poverty'), each = 6),
  year = rep(c(1970, 1980, 1990, 2000, 2010, 2012), 2),  
  concentration = c(.125, .12, .14, .13, .145, .146, .068, .09, .125, .119, .13, .135)
)

t1 <- textGrob(expression("Concentration of " * phantom(bold("affluence")) * "and" * phantom(bold("poverty")) * " nationwide"),
               x = 0.5, y = 1.1, gp = gpar(col = "purple"))

t2 <- textGrob(expression(phantom("Concentration of ") * bold("affluence") * phantom(" and poverty nationwide")),
               x = 0.5, y = 1.1, gp = gpar(col = "#EEB422"))

t3 <- textGrob(expression(phantom("Concentration of affluence and ") * bold("poverty") * phantom(" nationwide")),
               x = 0.5, y = 1.1, gp = gpar(col = "#238E68"))

# plot and add grobs with annotation_custom
ggplot(data, aes(year, concentration, color = group)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  annotation_custom(grobTree(t1, t2, t3)) +
  scale_y_continuous(limits = c(0, 0.15)) +
  scale_color_manual(values = c("#EEB422", "#238E68")) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL)+
  theme_minimal() +
  theme(legend.position = 'none',
        # add some extra margin on top
        plot.margin = unit(c(4, 1, 1, 1), "lines"))








