library(data.table)
library(RCurl)
library(stringr)
library(Sigfried)
library(reshape2)
library(rapportools)

#Read in gene annotation to map ENSEMBL gene IDs to Gene Symbols via GDSx WS: 
gencode_annotation <- read.csv('https://gdsx.prd.nibr.novartis.net:4000/getAnnotationData?annotation_name=gencode.gene&column_names=GENE_SYMBOL&media_type=text%2Fcsv&streaming=false')
names(gencode_annotation) <- c('ENSEMBL_ID', "GENE_SYMBOL")

sample_annotation <- readRDS('iobp_subset_noNormals_annotation.RDS')

getGdsxWsCall <- function(geneDataTable, observationTables) {
  base_url <- "https://gdsx.prd.nibr.novartis.net:4000/getData?"
  if (length(geneDataTable$ENSEMBL_ID) == 1) {
    ensembl=geneDataTable$ENSEMBL_ID
    gene_url=paste0("gene_ids=",ensembl)
    operator_head <- "&operator=ALL&head=false&head_size=5&media_type=text%2Fcsv"
    observation_url <- c()
    for (i in 1:length(observationTables$dataset_names)) {
        observation=observationTables[i]$gdsx_input_names
        observation_id=paste0("&observation_names=",observation)
        observation_url=paste0(observation_id, observation_url)
      }
    gdsx_ws_call <- paste0(base_url, gene_url, operator_head, observation_url)
  } else {
    genes_url <- c()
    geneDataTable <- as.data.table(geneDataTable)
    for (i in 1:length(geneDataTable$ENSEMBL_ID)) {
      ensembl=geneDataTable[i]$ENSEMBL_ID
      gene_ids=paste0("&gene_ids=",ensembl)
      genes_url = paste0(gene_ids, genes_url)
    }
    genes_url <- substring(genes_url,2)
    operator_head <- "&operator=ANY&head=false&head_size=5&media_type=text%2Fcsv"
    gdsx_ws_call <- paste0(base_url, genes_url, operator_head)
    observation_url <- c()
    for (i in 1:length(observationTables$dataset_names)) {
    observation=observationTables[i]$gdsx_input_names
    observation_id=paste0("&observation_names=",observation)
    observation_url=paste0(observation_id, observation_url)
    }
    gdsx_ws_call <- paste0(base_url, genes_url, operator_head, observation_url)
  }
  return(gdsx_ws_call)
}
 

fpkmToTpm <- function(fpkm)
 {
     exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

normalizeDatasets <- function(user_data, genelists) {
  tpm_datasets <- c('ptx.rnaseq.gene.tpm_5.1.0', 'ccle.rnaseq.gene.tpm_1.0.0', 'MMRF_CoMMpass_TPM_rnaseq_1.0.0')
  fpkm_datasets <- subset(user_data, !(OBSERVATION_NAME %in% tpm_datasets))
  print(unique(fpkm_datasets$OBSERVATION_NAME))
  if (any(unique(user_data$OBSERVATION_NAME) %in% tpm_datasets)) {
    tpm_datasets <- subset(user_data, OBSERVATION_NAME %in% tpm_datasets)
    tpm_datasets$VALUE <- as.numeric(tpm_datasets$VALUE)
    tpm_datasets <- dcast(tpm_datasets,SAMPLE_ID ~ GENE_ID, value.var="VALUE", fun.aggregate = mean)
    rownames(tpm_datasets) <- tpm_datasets$SAMPLE_ID
    tpm_datasets <- as.data.frame(data.matrix(subset(tpm_datasets, select=-(SAMPLE_ID))))
  }
  else {
    tpm_datasets <- c()
  }
  if (any(unique(user_data$OBSERVATION_NAME) %in% fpkm_datasets$OBSERVATION_NAME)) {
    fpkm_datasets$VALUE <- as.numeric(fpkm_datasets$VALUE)
    dataset <- dcast(fpkm_datasets, GENE_ID ~ SAMPLE_ID, value.var="VALUE", fun.aggregate = mean)
    rownames(dataset) <- dataset$GENE_ID
    dataset <- dataset[ , -which(names(dataset) %in% c("GENE_ID"))]
    dataset_mat <- data.matrix(dataset)
    dataset_tpm <- apply(dataset_mat, 1, function(x) fpkmToTpm(x))
    tpm_df <- as.data.frame(dataset_tpm)
    tpm_converted_data <- tpm_df
  }
  else {
    tpm_converted_data <- c()
  }
  if (!is.null(tpm_datasets) && !is.null(tpm_converted_data)) {
    all_tpm_data <- rbind(tpm_datasets, tpm_converted_data)
  } else if (!is.null(tpm_datasets) && is.null(tpm_converted_data)) {
    all_tpm_data <- tpm_datasets
  } else { 
    all_tpm_data <- tpm_converted_data
  }
  all_tpm_mat <- data.matrix(all_tpm_data)
  all_tpm_data_logged <- as.data.frame(apply(all_tpm_mat, 2, function(x) log2(x+1)))
  geneLists <- subset(genelists, genelists$ENSEMBL_ID %in% colnames(all_tpm_data_logged))
  names(all_tpm_data_logged) <- geneLists$GENE_SYMBOL
  all_tpm_data_logged$SampleID <- rownames(all_tpm_data_logged)
  return(all_tpm_data_logged)
}
#add_signatures
addCustomSignature <- function(normalized_dataset, custGenes, customSignatureName) {
  sig_mat <- data.matrix(normalized_dataset[ names(normalized_dataset)[names(normalized_dataset) %in% custGenes] ])
  cust_sig_calc <- apply(sig_mat, 1, function(x) mean((x)))
  cust_sig_calc <- as.data.frame(cust_sig_calc)
  names(cust_sig_calc) <- customSignatureName
  cust_sig_calc$SampleID <- rownames(cust_sig_calc)
  return(cust_sig_calc)
}

getFinalDataset <- function(normalized_data, cust_sig_calc) {
 
  if (is.empty(cust_sig_calc)) {
    final_dataset <- merge(sample_annotation, normalized_data, by=c('SampleID'))
  } else {
    normalized_data_custSig<- merge(normalized_data, cust_sig_calc)
    final_dataset <- merge(sample_annotation, normalized_data_custSig, by=c('SampleID'))
    return(final_dataset)
  }
}

