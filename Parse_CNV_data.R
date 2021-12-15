

## load libraries
library(maftools)
library(tidyverse)
library('biomaRt')
library(rjson) 

### ### ### START ### ### ### 

###### Part 1. Read CNV files

#####
## A) Colorectal
COAD <- read.csv("./Data/Colorectal/CNV/gdc_download_20200922_214647.478142/7f01e47a-2f4c-4b91-8db6-e1b6b5595390/COAD.focal_score_by_genes.txt", sep="\t")   
READ <- read.csv("./Data/Colorectal/CNV/gdc_download_20200922_214647.478142/8dcb2bfd-0e8a-402d-941c-def06f2e8a96/READ.focal_score_by_genes.txt", sep="\t")   
CNVdata <- merge(COAD, READ, by = c('Gene.Symbol', "Gene.ID", "Cytoband"))
colnames(CNVdata) # 679 -3 = 676 samples

#take sample names
CNVfilenames <-  colnames(CNVdata) 
CNVfilenames <- CNVfilenames[-c(1,2,3)]

##### Part 2. Map gene ids to ensembl ids

# Make a Biomart query

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- CNVdata$Gene.Symbol
G_list <- getBM(attributes= c("ensembl_gene_id", "hgnc_symbol"),values=genes,mart= mart)  
# I used gene_id again. 


# The ens ID are a bit different from the ensembl's. To make them the same we keep the first 15 chars on the names.
ens_id <- CNVdata$Gene.Symbol
ens_id <- sapply(ens_id, function(x){substring(x, 1, 15)}) # take 15 first elements of each string/name
CNVdata$Gene.Symbol <- ens_id # change the gene names in our dataframe

CNVdata <- merge(CNVdata,G_list,by.x = "Gene.Symbol", by.y="ensembl_gene_id") #gerge data with gene symbols

colnames(G_list) <- c("Gene.Symbol", "hgnc_symbol")
a <- merge(CNVdata,G_list,by = "Gene.Symbol") #gerge data with gene symbols


CNVdata$Gene.Symbol <- NULL # Delete ensembe id col
CNVdata <- as.tibble(CNVdata)%>% dplyr::select(hgnc_symbol, everything()) ##### 
CNVdata <- CNVdata %>% mutate_all(na_if,"") # put NA to enseble ids that dont map to our genes #19453 rows
CNVdata <- na.omit(CNVdata) # delete rows that dont map to genes 19061 rows

##### Part 3. Select only genes needed (Wnt/b-catenin related genes)

#Load the list of gene we will use
cangenes <- read.csv("APC_reporter_gene_candidates.txt", header = FALSE)
CNVdata <- CNVdata[CNVdata$hgnc_symbol %in% cangenes$V1,] # take only the genes we want 350 obs
names(CNVdata)[1] <- "Genes"  # rename hgnc symbol to Genes

# transpose the table so we have genes as cols and samples as rows
CNVdata <- as.data.frame(t(CNVdata))
names(CNVdata) <- lapply(CNVdata[1, ], as.character) #put correct col names
CNVdata <- CNVdata[-1,] 
CNVdata <- tibble::rownames_to_column(CNVdata, "Samples") # make the row names into first col and name it genes

#remove and keep first 2 rows
firtstrows <- CNVdata[1:2,] ### Question: Do I have to place them back in? Its Gene.ID and Cytoband
CNVdata <- CNVdata[-c(1:2),] 

#Replace file names xx.xx to xx_xx to match 
CNVdata$Samples <- gsub('\\.', '-', CNVfilenames)


##### Part 4. Find sample/ barcode match

## Load CNV metadata and take the match of samples with barcodes
metadata <- fromJSON(file = "./Data/Colorectal/CNV/CNV_metadata.json") 

####
## A) Colorectal, Esophagogastric
COADmetadata <- metadata[[2]][["associated_entities"]]
COADname <- sapply(COADmetadata, function(x){x$entity_id[[1]]})
COADbc <- sapply(COADmetadata, function(x){x$entity_submitter_id[[1]]})

READmetadata <- metadata[[1]][["associated_entities"]]
READname <- sapply(READmetadata, function(x){x$entity_id[[1]]})
READbc <- sapply(READmetadata, function(x){x$entity_submitter_id[[1]]})

# Sample-bc match dataframe
CNVmatch <- data.frame(c(COADname, READname),c(COADbc, READbc))
names(CNVmatch) <- c("Samples", "barcode")

## B) Endometrial
UCECmetadata <- metadata[[1]][["associated_entities"]]
name <- sapply(UCECmetadata, function(x){x$entity_id[[1]]})
bc <- sapply(UCECmetadata, function(x){x$entity_submitter_id[[1]]})

# Sample-bc match dataframe
CNVmatch <- data.frame(name,bc)
names(CNVmatch) <- c("Samples", "barcode")
######

## Some samples have an "X" in from of their names so they dont match
# We have to remove that X
interdf <- as.data.frame(CNVmatch$Samples)
interdf$CNVdataSamples <- CNVdata$Samples

# find which samples statr with the "X"
thex<- substr(CNVdata$Samples, 0, 1) == "X" #which samples have the x in the front
interdf$CNVdataSamples[thex] <- sapply(interdf$CNVdataSamples[thex], function(x){sub('.', '', x)}) # find the samples and delete first letter
CNVdata$Samples <- interdf$CNVdataSamples
## Match the data
#use the interdf
int <- intersect(CNVmatch$Samples, CNVdata$Samples) # all 676 samples match

#merge the data frames and put corret labels
jointdataset <- merge(CNVdata, CNVmatch, by = 'Samples') #merge the data frames

columnNumber <- which(colnames(jointdataset)=="barcode")
jointdataset <- jointdataset[,c(columnNumber,1:ncol(jointdataset)-1)]
jointdataset$Samples <- NULL

# write table in file
write.csv(jointdataset,"./Data/Colorectal/CNV/CRC_CNV_data_prematched", row.names = FALSE)

#### ### ### END ### ### ###

