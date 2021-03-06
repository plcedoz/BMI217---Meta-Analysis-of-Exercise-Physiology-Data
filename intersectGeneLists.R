#This is a function to intersect our computed gene lists with external gene lists of interest: ELITE, HCM, LQTS and other
#cardiovascular diseases

dataPath <- "Data/"

#paths to gene lists
#all gene lists are in list format with one gene per line
#cardiome: genes relevant for cardiac phenotypes
cardiomePath <- paste(dataPath, 'cardiome_genes_list.txt.updated', sep = '')
#sarcomere: sarcomeric related genes
sarcomerePath <- paste(dataPath, "sarcomericGenes.txt", sep = '')
#clinvar: genes associated with disease pathogenicity
clinvarPath <- paste(dataPath, "clinvar_gene.txt", sep = '')
#elite lof : loss of function genes identified as significant in the elite study
elitePath <- paste(dataPath, "elite_lof_genes.txt", sep = '')

#GENE dx data-->cardiac disease associated genes
#arvc (arrythmogenic right ventricular cardiomyopathy)
arvcPath <- paste(dataPath, "genedx/arvc.txt", sep = '')
#brugada (syndrome associated with sudden death related to cardiac phenotypes)
brugadaPath <- paste(dataPath, "genedx/brugada.txt", sep = '')
#cardioFacioCutaneous (a disorder assocaited with weird things in the heart and skin)
cardioFacioCutaneousPath <- paste(dataPath, "genedx/cardio-facio-cutaneous.txt", sep = '')
#cardiomyopathy genes
cardioMyopathyPath <- paste(dataPath, "genedx/comprehensive_cardiomyopathy.txt", sep = '')
#cpvt: catecholaminergic polymorphic ventricular tachycardia
cpvtPath <- paste(dataPath, "genedx/cpvt.txt", sep = '')
#dcm: dialted cardiomyopathy
dcmPath <- paste(dataPath, "genedx/dcm.txt", sep = '')
#hcm: hypertrophies cardiomyopathy
hcmPath <- paste(dataPath, "genedx/hcm.txt", sep = '')
#lds: loeys-dietz syndrome (cardiac tissue syndrome)
ldsPath <- paste(dataPath, "genedx/lds.txt", sep = '')
#lqts: long qt syndrome
lqtsPath <- paste(dataPath, "genedx/lqts.txt", sep = '')
#lvnc: left non-ventricular contract cardiomyopathy
lvncPath <- paste(dataPath, "genedx/lvnc.txt", sep = '')
#noonan syndrome 
noonanPath <- paste(dataPath, "genedx/noonan.txt", sep = '')

#read in our gene lists (note they are in heterogeneous formats so they need to be read in different formats

#bad heart phenotype related genes
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
lvnc[2] = "ACTC"
lvnc[27] = "LDB3"
noonan <- read.csv(noonanPath, sep = ',', header = FALSE)
#again apply the correction for string matching
noonan[2] = "ACTC"
noonan[27] = "LDB3"

#CLINVar 
clinvar <- read.csv(clinvarPath, sep = '\t', header = TRUE)
clinvar <- clinvar$gene
clinvar <- as.character(clinvar)

#ELITE LOF
elite <- read.csv(elitePath, sep = '\t', header = FALSE)


athlete_poster <- read.csv("Data/AlthleteGenes.txt", sep = ',', header = FALSE)


#create a matrix of genelists
geneLists <- list(arvc, brugada, cardiomyopathy, cpvt, dcm, hcm, lds, lqts, lvnc, noonan, clinvar, elite, athlete_poster)    
#create the names for this genelist (make sure they match! EVERYTHING IS RUINED IF THEY DONT MATCH)
names(geneLists) <- c("arvc", "brugada", "cardiomyopathy", "cpvt", "dcm", "hcm", "lds", "lqts", "lvnc", "noonan", "clinvar", "elite", "athlete_poster")

#given a primary gene list, intersects it with each list on the list of lists and writes a file with the name associated with the file enumerated in the list of lists
intersect_and_write_lists <- function(primaryGeneList, listOfLists){
 for(i in 1:length(listOfLists)){
   name <- names(listOfLists)[i]
   curList <- as.matrix(get(name, listOfLists))
   #we pass primary gene list second because its gene names may be gross
   l <- intersect(curList, primaryGeneList)
   #ALERT! you should change this path to where you want the files to be written
   writePath <- paste("output_data/intersection_with_", name, ".txt", sep = "")
   write(as.matrix(l), writePath)
 } 
}

intersect_and_write_lists(final_gene_list, geneLists)




