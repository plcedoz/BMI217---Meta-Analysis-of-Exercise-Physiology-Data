#paths to gene lists
#all gene lists are in list format with one gene per line
#cardiome: genes relevant for cardiac phenotypes
cardiomePath <- "/Users/noahfriedman/Desktop/BMI217Project/data/cardiome/cardiome_genes_list.txt.updated"
#sarcomere: sarcomeric related genes
sarcomerePath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/sarcomericGenes.txt"

#GENE dx data-->cardiac disease associated genes
#arvc (arrythmogenic right ventricular cardiomyopathy)
arvcPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/arvc.txt"
#brugada (syndrome associated with sudden death related to cardiac phenotypes)
brugadaPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/brugada.txt"
#cardioFacioCutaneous (a disorder assocaited with weird things in the heart and skin)
cardioFacioCutaneousPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/cardio-facio-cutaneous.txt"
#cardiomyopathy genes
cardioMyopathyPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/comprehensive_cardiomyopathy.txt"
#cpvt: catecholaminergic polymorphic ventricular tachycardia
cpvtPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/cpvt.txt"
#dcm: dialted cardiomyopathy
dcmPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/dcm.txt"
#hcm: hypertrophies cardiomyopathy
hcmPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/hcm.txt"
#lds: loeys-dietz syndrome (cardiac tissue syndrome)
ldsPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/lds.txt"
#lqts: long qt syndrome
lqtsPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/lqts.txt"
#lvnc: left non-ventricular contract cardiomyopathy
lvncPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/lvnc.txt"
#noonan syndrome 
noonanPath <- "/Users/noahfriedman/Desktop/BMI217Project/data/genesets/genedx/noonan.txt"

#read in our gene lists (note they are in heterogeneous formats so they need to be read in different formats
arvc <- read.csv(arvcPath, sep = ',', header = FALSE)
brugada <- read.csv(brugadaPath, sep = ',', header = FALSE)
cardiomyopathy <- as.data.frame(t(read.csv(cardioMyopathyPath, sep = '\n', header = FALSE)))
cpvt <- read.csv(cpvtPath, sep = ',', header = FALSE)
dcm <- read.csv(dcmPath, sep = ',', header = FALSE)
#we fix the second entry so it says ACTC instead of ACTC(ACTC1)
dcm[2] = "ACTC"
hcm <- as.data.frame(t(read.csv(hcmPath, sep = '\n', header = FALSE)))
#again we fix two weird entries
hcm[2] = "ACTC"
hcm[28] = "LDB3"
lds <- read.csv(ldsPath, sep = ',', header = FALSE)
lqts <- as.data.frame(t(read.csv(lqtsPath, sep = '\n', header = FALSE)))
lvnc <- read.csv(lvncPath, sep = ',', header = FALSE)
#again apply the correction for string matching
#TEST
lvnc[2] = "BUTTFACE"
lvnc[27] = "LDB3"
noonan <- read.csv(noonanPath, sep = ',', header = FALSE)
#again apply the correction for string matching
noonan[2] = "ACTC"
noonan[27] = "LDB3"

#create a matrix of genelists
geneLists <- list(arvc, brugada, cardiomyopathy, cpvt, dcm, hcm, lds, lqts, lvnc, noonan)    
#create the names for this genelist (make sure they match! EVERYTHING IS RUINED IF THEY DONT MATCH)
names(geneLists) <- c("arvc", "brugada", "cardiomyopathy", "cpvt", "dcm", "hcm", "lds", "lqts", "lvnc", "noonan")

#gets the intersection of two gene lists, gl1, and gl2
#because it uses grepl, you must put the shorter gene name first and the longer gene name last
get_intersection <- function(gl1, gl2){
  tfVec <- rep(0, length(gl1))
  for(i in 1: length(gl1)){
    include <- FALSE
    gl1Name <- gl1[i]
    v <- grepl(gl1Name, gl2)
    if(TRUE %in% v){
      include <- TRUE
    }
    else{
      print("ochen ploxa")
      print(gl1Name)
    }
    tfVec[i] = include
  }
  tfVec <- as.logical(tfVec)
  #we return the original list gl1, but subsetted only for intersections
  #there's a bit of rigamarole with turning things into an out of dataframes here
  return(as.data.frame(gl1[tfVec]))
}

#given a primary gene list, interests it with each list on the list of lists and writes a file with the name associated with the file enumerated in the list of lists
intersect_and_write_lists <- function(primaryGeneList, listOfLists){
 for(i in 1:length(listOfLists)){
   name <- names(listOfLists)[i]
   curList <- get(name, listOfLists)
   #we pass primary gene list second because its gene names may be gross
   l <- get_intersection(curList, primaryGeneList)
   writePath <- paste("/Users/noahfriedman/Desktop/BMI217Project/INTERSECTIONS_WITH", name, ".txt", sep = "")
   write(as.matrix(l), writePath)
 } 
}

intersect_and_write_lists(geneList2$Gene.symbol, geneLists)



