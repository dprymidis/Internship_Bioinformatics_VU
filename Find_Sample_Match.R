
### READ ALL DATA IN 1 SCRIPT
# Finds MSS samples with SNV, CNV and gene expression data
# Selects these samples in SNV, CNV, gene expression files and generate new ones
# CNV and Gne_Expression Done
# SNV failed


## load libraries
library(maftools)
library(tidyverse)
library(rjson) 


### ### ### START ### ### ###

#### Part 1. read files 

## read SNV files
## A) COAD/READ, EASCA/STAD
coadSNVmaf <- read.maf("./Data/Colorectal/SNV/gdc_download_20200926_134053.784880/03652df4-6090-4f5a-a2ff-ee28a37f9301/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.gz")
readSNVmaf <- read.maf("./Data/Colorectal/SNV/gdc_download_20200926_134053.784880/faa5f62a-2731-4867-a264-0e85b7074e87/TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf.gz")
SNVdata <- merge_mafs(c(coadSNVmaf, readSNVmaf))

## B) UCEC,
  SNVdata <- read.maf("./Data/Endometrial/SNV/gdc_download_20200925_180842.663826/d3fa70be-520a-420e-bb6d-651aeee5cb50/TCGA.UCEC.mutect.d3fa70be-520a-420e-bb6d-651aeee5cb50.DR-10.0.somatic.maf.gz")
## C) Egophagogastric
## D) Melanoma
SNVdata <- read.maf("./Data/Melanoma/SNV/gdc_download_20200925_082238.004247/4b7a5729-b83e-4837-9b61-a6002dce1c0a/TCGA.SKCM.mutect.4b7a5729-b83e-4837-9b61-a6002dce1c0a.DR-10.0.somatic.maf.gz")
## E) Non-small Cell Lung
LUADSNVmaf <- read.maf("./Data/Non_Small_Cell_Lung/SNV/gdc_download_20200925_204926.413520/0458c57f-316c-4a7c-9294-ccd11c97c2f9/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf.gz")
LUSCSNVmaf <- read.maf("./Data/Non_Small_Cell_Lung/SNV/gdc_download_20200928_124012.219988/95258183-63ea-4c97-ae29-1bae9ed06334/TCGA.LUSC.mutect.95258183-63ea-4c97-ae29-1bae9ed06334.DR-10.0.somatic.maf.gz")
SNVdata <- merge_mafs(c(LUADSNVmaf, LUSCSNVmaf))
######

# read CNV files
## Load CNV data and take barcodes
CNVdata <- as.data.frame(read.csv("./Data/Colorectal/CNV/CRC_CNV_data_prematched"))
names(CNVdata)[1] <- "Tumor_Sample_Barcode"

## read gene expresion files 
# Counts
Genedata <- as.data.frame(read.csv("./Data/Colorectal/Gene_expression/CRC_Gene_expression_all_Counts"))
names(Genedata)[1] <- "Tumor_Sample_Barcode"  # rename first col

## load MSS information, only for colorectal and esophagogastric
obsMSS <- as.data.frame(read.csv("./CRC_MSS_samples"))
obsMSS$X <- c()
names(obsMSS)[1] <- "Tumor_Sample_Barcode"

###### Part 2. Find samples with all: MSS, CNV, SNC and gene data  (using barcodes)

# Check num of observations
obsSNV <- SNVdata@variants.per.sample[,1]
obsCNV <- as.data.frame(CNVdata$Tumor_Sample_Barcode)
names(obsCNV)[1] <- "Tumor_Sample_Barcode"
obsGene <- as.data.frame(Genedata$Tumor_Sample_Barcode)
names(obsGene)[1] <- "Tumor_Sample_Barcode"


# make list of barcodes
# take 15 first elements of each string/name because the first 16 letters define a sample per parient
Sbc <-as.vector(sapply(obsSNV$Tumor_Sample_Barcode, function(x){substring(x, 1, 15)}))  ### 535
Cbc <-as.vector(sapply(obsCNV$Tumor_Sample_Barcode, function(x){substring(x, 1, 15)})) # 536
Gbc <-as.vector(sapply(obsGene$Tumor_Sample_Barcode, function(x){substring(x, 1, 15)})) # 604
Mbc <- as.vector(sapply(obsMSS$Tumor_Sample_Barcode, function(x){substring(x, 1, 15)})) # 444


# check for doublecats
duplicated(Sbc) 
duplicated(Cbc) 
Cbc <- Cbc[!duplicated(Cbc)] #### 590
duplicated(Gbc) 
Gbc <- Gbc[!duplicated(Gbc)] #### 590
duplicated(Mbc) 


### check the SNV and CNV

CSbc <- intersect(Cbc,Sbc) # 522 samples

## chech with Gene data
GCSbc <- intersect(CSbc,Gbc) # 449 samples

# check with MSS data
MGCSbc <- intersect(GCSbc,Mbc) # 386 samples

Final_Samples <- as.data.frame(MGCSbc)
names(Final_Samples)[1] <- "Tumor_Sample_Barcode"

write.csv(Final_Samples,"./CRC_Matched_Samples", row.names =TRUE)

### Part 3. Select the data from samples and write new data files

## For gene expression
# Change barcodes in data
Genedata$Tumor_Sample_Barcode <- sapply(Genedata$Tumor_Sample_Barcode, function(x){substring(x, 1, 15)})

# Remove duplicates
Genedata <- Genedata[!duplicated(Genedata$Tumor_Sample_Barcode),] # from 604 goes to 590
Final_Gene_data <- Genedata %>% filter(Tumor_Sample_Barcode  %in% Final_Samples$Tumor_Sample_Barcode)

# Write the gene data in a file
write.csv(Final_Gene_data,"./Data/Colorectal/Gene_expression/CRC_Gene_expression_Counts", row.names = TRUE)


## For CNV data
CNVdata$Tumor_Sample_Barcode <- sapply(CNVdata$Tumor_Sample_Barcode, function(x){substring(x, 1, 15)})

# Remove duplicates
CNVdata <- CNVdata[!duplicated(CNVdata$Tumor_Sample_Barcode),] # from 676 goes to 620
Final_CNV_data <- CNVdata %>% filter(Tumor_Sample_Barcode  %in% Final_Samples$Tumor_Sample_Barcode)

# Write the gene data in a file
write.csv(Final_Gene_data,"./Data/Colorectal/CNV/CRC_CNV_data_matched", row.names = TRUE)

## For SNV data
#I cant yet. I use the matched samples file to parse the samples from the maf files
########## END ##########

