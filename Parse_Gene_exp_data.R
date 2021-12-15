
## In order to run the script for ie. Counts and CRC
# First extract the data 
# Then execute this command in the terminal:
# find ~/Documents/Master/internship/APC/APC_Project/Data/Colorectal/Gene_expression/Counts/gdc_download_20200907_123430.425003/ -type f -print0 | xargs -0 mv -t ~/Documents/Master/internship/APC/APC_Project/Data/Colorectal/Gene_expression/Counts/
# This command takes all the files from the subfolders and places them in the main folder


### Load libraries
library(tidyverse)
library('biomaRt')
library(rjson) 

### ### ### Start ### ### ###

#### part 1. Load gene expression data 

# List all file names
namelist = list.files(path = "./Data/Colorectal/Gene_expression/Counts/", pattern = "*counts.gz") 

# This doesnt work but I dont know why and I had to change working dir
# datalist = lapply(namelist, function(x)read.table(x))  #make a table of contents of files (nested list)
setwd("~/Documents/Master/internship/APC/APC_Project/Data/Colorectal/Gene_expression/Counts/")
datalist = lapply(namelist, function(x)read.table(x))  #make a table of contents of files (nested list)
setwd("~/Documents/Master/internship/APC/APC_Project") # return 

# When we load the data the first col (with the segment names) is repeated. we only the the 1st on time and then all the 2nd cols
rnames <-datalist[[c(1,1)]] #take row names
cnames <- sapply(datalist, "[", 2) #take all data (2nd cols) without the row names

# Put all data into a data frame
genexpdata = do.call("cbind", cnames) 
row.names(genexpdata) <- rnames # Correct names in rows
colnames(genexpdata)  <- namelist

##### Part 2. Map gene ids to ensembl ids

# need to install biomaRT to map the ensembl IDs to gene names
# BiocManager::install("biomaRt")

genexpdata <- as.data.frame(genexpdata) # make genexp a data frame for better use
genexpdata <- tibble::rownames_to_column(genexpdata, "genes") # make the row names into first col and name it genes

# Make a Biomart query
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- genexpdata$genes
G_list <- getBM(attributes= c("ensembl_gene_id", "hgnc_symbol"),values=genes,mart= mart)   #ensemble_gene_id or ensemble_transcript_id? 
# I get a warning
# I used gene_id because it had number of observations close to our list. 
# When transcript id was used the obs were much more

# Our ens ID are a bit different from the ensembl's. To make them the same we keep the first 15 chars on the names.
ens_id <- genexpdata$genes
ensem_ids <- sapply(ens_id, function(x){substring(x, 1, 15)}) # take 15 first elements of each string/name
genexpdata$genes <- ensem_ids # change the gene names in our dataframe
maxgenedata <- merge(genexpdata,G_list,by.x = "genes", by.y="ensembl_gene_id") #put everything into the dataframe

maxgenedata$genes <- NULL # Delete ensembe id col
maxgenedata <- as.tibble(maxgenedata)%>% dplyr::select(hgnc_symbol, everything()) # Places hgnc_symbol 1st col
# warning again
maxgenedata <- maxgenedata %>% mutate_all(na_if,"") # put NA to enseble ids that dont map to our genes #56501 rows
genedata <- na.omit(maxgenedata) # delete rows that dont map to genes 37808 rows
Genedata <- genedata

# transpose the table so we have genes as cols and samples as rows
Genedata <- as.data.frame(t(Genedata))
names(Genedata) <- lapply(Genedata[1, ], as.character) #put correct col names
Genedata <- Genedata[-1,] 
Genedata <- tibble::rownames_to_column(Genedata, "file") # make the row names into first col and name it genes

##### Part 4. #####
# load metadata in order to connect file names to sample barcodes
# For Counts, FPKM and FPKM-UQ adn for each cohort we use different metadata
metadata <- fromJSON(file = "Data/Colorectal/Gene_expression/Gene_expression_metadata-_counts.json") # gene expression metadata

allnames <- sapply(metadata, function(x){x$file_name}) #take all file
allbcs <- sapply(metadata, function(x){x$associated_entities[[1]]$entity_submitter_id}) # take all sampe barcodes # old
matchnames <- as.data.frame(allnames, allbcs)
matchnames <- do.call(rbind, Map(data.frame, A=allnames, B=allbcs))
colnames(matchnames) <- c("file", "barcode")

#merge the data frames and put corret labels
jointdataset <- merge(Genedata, matchnames, by = 'file') #merge the data frames
jointdataset <- tibble::column_to_rownames(jointdataset, var = "barcode") #put barcodes as row names
jointdataset$file <- c() # delete file names

# we have the gene data
head(jointdataset)

write.csv(jointdataset,"./Data/Colorectal/Gene_expression/CRC_Gene_expression_all_Counts", row.names = TRUE)
### ### ### END ### ### ###


#### Check sample tissue sourse
matchnames$barcode
# find which samples statr with the "X"
sort(unique(substr(matchnames$barcode, 6, 7))) #find all unique tissue sourses in order
 
### ### ### END ### ### ###

## A) Colorecta
## B) Endometrial
## C) Esophagogastric
## D) Melanoma
## E)Non_Small_Cell_Lung

