
## load libraries
library(maftools)
library(tidyverse)
library(BiocManager)
library(TCGAbiolinks)


# How to install packages
# install.packages("BiocManager")
# BiocManager::install("TCGAbiolinks")

### ### ### Start ### ### ###

### Part 1. Load the data

#####
## A) Colorectal
COAD <- read.maf("./Data/Colorectal/SNV/gdc_download_20200926_134053.784880/03652df4-6090-4f5a-a2ff-ee28a37f9301/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.gz")   
READ <- read.maf("./Data/Colorectal/SNV/gdc_download_20200926_134053.784880/faa5f62a-2731-4867-a264-0e85b7074e87/TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf.gz")   
SNVdata <- merge_mafs(c(COAD, READ))


### ### Part 2. Plot data

# Variat number and sample ID 
samplevar <- SNVdata@variants.per.sample 
samplevar$log10var = log10(samplevar$Variants)   
samplevar$plotpoint = 10    
set.seed(9)

# creates random numbers between 10 and 14 for plot
randomnum <- runif(536, 10, 14) 

# assign random numbers to samples randomly
samplevar$plotpoint <- randomnum 

#### plot samples and variants
p1 <- ggplot(data=samplevar, mapping=aes(x=plotpoint, y=log10var)) +
  geom_point(size=2) + 
  ggtitle("Tumor Mutational Burden") +
  labs(main = "Tumor Mutational Burden" , x = "Samples", y = "Log10 variant number") +
  xlim(9, 15)
print(p1)

###remove obvious outliers 
samplevar <- filter(samplevar, !Tumor_Sample_Barcode== "TCGA-AH-6547-01A-11D-1826-10")


### Plot Sample mutational burden distribution
p2 <- ggplot(data=samplevar, mapping=aes(x=plotpoint, y=log10var)) +   # nice plot
  geom_point(size=1) + 
  ggtitle("Colorectal cancer TMB") +
  labs(main = "Tumor Mutational Burden" , x = "Samples", y = "Log10 variant number") +  
 # geom_hline(yintercept = 2.51, color = "red", size = 1) +
  xlim(9, 15)
print(p2)


# Save the plot in a png file
#####
png("./CRC_TMB_Plot.png")
print(p2)
dev.off()


#### Part 3. Add MS status
# get MS status from GDC legacy data

###installation of packages
#install.packages("BiocManager")
#BiocManager::install("BioinformaticsFMRP/TCGAbiolinks") 

## COAD READ, UCEC, ESCA STAD, cant find SKCM, LUAD and LUSC
query <- GDCquery(project = c("TCGA-COAD", "TCGA-READ"),  
                  data.category = "Other",
                  legacy = TRUE,
                  access = "open",
                  data.type = "Auxiliary test")
GDCdownload(query)
msi_results <- GDCprepare_clinic(query, "msi")
colnames(msi_results)[3] <- "MS_Status"


# In the MSI_reults they use the 12 first letters of the tumor sample barcode 
x <- samplevar$Tumor_Sample_Barcode  # we want the 12 first letters in order to match
samplevar$bcr_patient_barcode <- substring(x, 1, 12)

#merge the two data frames te names with logvar and the msi
samplevar <- merge(samplevar, msi_results, by = 'bcr_patient_barcode')  

# Plot log10 sample variance with MS status
p4 <- ggplot(data=samplevar, mapping=aes(x=plotpoint, y=log10var)) +   ##### plot all of our samples
  geom_point(aes(colour = MS_Status), size = 1) + 
  ggtitle("MS Status in Colorectal cancer") +
  labs(main = "Tumor Mutational Burden" , x = "Samples", y = "Log10 variant number") +  
  geom_hline(yintercept = 2.5, color = "red", size = 1) +
  xlim(9, 15)
print(p4)

# Save the plot in a png file
png("./MS_Plot.png")
print(p4)
dev.off()


#### Part 4 Exclude MS high and POLE mut samples and write the barcodes of MSS samples in a text
# keep samples bellow 


samples_used <- (samplevar$log10var < 2.5)
final_samples <- samplevar[samples_used]  

write.csv(final_samples$Tumor_Sample_Barcode,"./CRC_MSS_samples", row.names = TRUE)

### ### ### ### ### END ### ### ### ### ###

