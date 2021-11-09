## Initialzie
library(GEDI) # Load the GEDI package
PATH_TO_DATA_FOLDERS <- "man/examples" # In github repository

# Read files with variables used in the functions
data_info <- read.table(paste(PATH_TO_DATA_FOLDERS,"/data_info.txt",sep = ""),
                        header = TRUE)
pheno <- read.table(paste(PATH_TO_DATA_FOLDERS,"/pheno.txt",sep = ""),
                    header = TRUE)

## 1. Read datasets with ReadGE
datasets <- ReadGE(dataFolders = data_info$dataFolders,
                   sources = data_info$sources,
                   path = PATH_TO_DATA_FOLDERS)

## 2. Integrate datasets
dat <- GEDI(datasets = datasets,
            attributes = data_info$attributes,
            species = "btaurus", # Bos taurus
            mapTo = "ensembl_gene_id",
            path = PATH_TO_DATA_FOLDERS)

## 3. Remove batch effect
cData <- BatchCorrection(Data = dat,
                         batch = pheno$batch,
                         status = pheno$status,
                         visualize = TRUE)

## 4. Verify data integration
res <- VerifyGEDI(X = cData,
                  y = pheno$status,
                  batch = pheno$batch,
                  model = "logistic")
