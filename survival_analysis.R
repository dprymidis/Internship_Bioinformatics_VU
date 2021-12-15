
library(readxl)
library(survival)
library(survminer)
library(tidyverse)

setwd("~/Documents/Master/internship/APC/APC_Project")
## Genomic_metrics is a dataframe with class labels as Class and the barcode as sample id. You can use your own class labels if you create such a dataframe
genome_metrics <- read.csv("./1_trunc_samples.filtered") 
colnames(genome_metrics) <- c("barcode","Class" )
genome_metrics$barcode <- sapply(genome_metrics$barcode, function(x){substring(x, 1, 12)})

# Remove CNA if you like
#CNV <- read.csv("./produced_data/APC_CNV") 
#colnames(CNV) <- c("cna", "Tumor_Sample_Barcode")
#CNV$Tumor_Sample_Barcode <- sapply(CNV$Tumor_Sample_Barcode, function(x){substring(x, 1, 12)})
#a<- CNV %>% filter(cna != 0)
#genome_metrics <- genome_metrics[!genome_metrics$barcode %in% a$Tumor_Sample_Barcode,]



##
# read clinical data, curated clinical data from Liu,J. et al. (2018) An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell, 173, 400-416.e11.
clinical_data = read_excel("mmc1.xlsx", sheet = "TCGA-CDR")

# process survival data, select overall survival, disease free interval, progression free interval and disease specific survival
survival_data = subset(clinical_data[which(clinical_data$bcr_patient_barcode %in% unique(genome_metrics$barcode)),], 
                       select = c(bcr_patient_barcode, OS.time, OS, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time))


# set numeric columns as numeric
survival_data = survival_data[!duplicated(survival_data$bcr_patient_barcode),]
survival_data$OS.time <- as.numeric(as.character(survival_data$OS.time))
survival_data$DSS.time <- as.numeric(as.character(survival_data$DSS.time))
survival_data$DFI.time <- as.numeric(as.character(survival_data$DFI.time))
survival_data$PFI.time <- as.numeric(as.character(survival_data$PFI.time))
survival_data[survival_data == "[Not Available]"] <- NA
survival_data$OS <- as.numeric(survival_data$OS)
survival_data$DSS <- as.numeric(survival_data$DSS)
survival_data$DFI <- as.numeric(survival_data$DFI)
survival_data$PFI <- as.numeric(survival_data$PFI)

# overall survival analysis
sfit = survfit(Surv(OS.time, OS)~1, data = survival_data)
plt = ggsurvplot(sfit)
plot.new()
print(plt, newpage = F)

# create TBL status
a<- genome_metrics[genome_metrics$barcode %in% survival_data$bcr_patient_barcode,]
Class<-as.factor(a$Class)

# PFI progresseion free, OS overall survival, DFI disease-free interval  , DSS Disease specific survival
# overall survival analysis with regard to the variable Class, this is a factor, you will get the survival of the different factor levels
sfit_Class = survfit(Surv(DSS.time, DSS) ~ Class, data = survival_data)
plt = ggsurvplot(sfit_Class, title = "Disease spesific survival")
plot.new()
print(plt, newpage = F)

# Cox regressional hazard model to calculate the hazard rate
cox_TBL <- coxph(Surv(OS.time, OS) ~ Class, data = survival_data)
