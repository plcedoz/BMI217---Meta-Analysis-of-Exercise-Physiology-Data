library(GEOquery)

get_data_gse <- function (gse_ref){
  gse <- getGEO(gse_ref, GSEMatrix=FALSE)
  meta <- Meta(gse)
  
  #List of samples
  gsms = names(GSMList(gse))
  
  #Extract charachteristics for all samples (gsms)
  gsm_description = c()
  for (gsm in gsms){
    charac = GSMList(gse)[gsm][[1]]@header$characteristics_ch1
    gsm_description = rbind(gsm_description, c(gsm, charac))
  }
  gsm_description = as.data.frame(gsm_description)
  
  #Verify that there is only one platform (GPL)
  GPLs = names(GPLList(gse))
  
  #Extract the probe/genes relationship from GPL
  #columns = GPLList(gse)[[1]]@dataTable@columns
  #datatable = GPLList(gse)[[1]]@dataTable@table
  probe_to_gene = GPLList(gse)[[1]]@dataTable@table$`Gene Symbol`
  
  #Create the matrix of samples vs genes
  probesets <- Table(GPLList(gse)[[1]])$ID
  data_matrix <- do.call("cbind", lapply(GSMList(gse), function(x) {
    tab <- Table(x)
    mymatch <- match(probesets, tab$ID_REF)
    return(tab$VALUE[mymatch])
  }))
  data_matrix <- apply(data_matrix, 2, function(x) {
    as.numeric(as.character(x))
  })
  
  #Do we need to log normalize?
  #Should we use quantile normalization?
  data_matrix <- log2(data_matrix)
  rownames(data_matrix) = probe_to_gene
  list(data_matrix = data_matrix, gsm_description = gsm_description)
}