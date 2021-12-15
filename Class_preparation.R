
### Prepare data for classification
## load libraries
library(maftools)
library(tidyverse)

### ### ### START ### ### ###

#### Part 1. Read files and take the samples we want

####
setwd("~/Documents/Master/internship/APC/APC_Project")
coadsnv <- read.table("./Data/Colorectal/SNV/gdc_download_20200926_134053.784880/03652df4-6090-4f5a-a2ff-ee28a37f9301/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.gz", sep = "\t" , header = TRUE, quote = "", stringsAsFactors = FALSE)
readsnv <- read.table("./Data/Colorectal/SNV/gdc_download_20200926_134053.784880/faa5f62a-2731-4867-a264-0e85b7074e87/TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf.gz", sep = "\t" , header = TRUE, quote = "", stringsAsFactors = FALSE)
SNVdata<- rbind(coadsnv,readsnv)

## Load the sample sheet of the according cohort
matchedsamples <- read.csv("./CRC_Matched_Samples")
matchedsamples[1] <- NULL
matchedsamples <- as.character(matchedsamples$Tumor_Sample_Barcode)

## change names in data to be the first 15 to make the selection 
SNVdata$Tumor_Sample_Barcode <- sapply(SNVdata$Tumor_Sample_Barcode, function(x){substring(x, 1, 15)})
SNVpost <- SNVdata[SNVdata$Tumor_Sample_Barcode %in% matchedsamples,] # take only MSS samples that have all the data

## Take all APC mutations data 
APCdata <- SNVpost %>% filter(Hugo_Symbol == "APC")
mutpos <- APCdata$Protein_position # take the positions
mutpos <- substr(mutpos, start = 1, stop = 4) %>% str_remove ("/") %>% str_remove ("-") # only the starting position
mutposdf <- as.data.frame(as.numeric(mutpos)) # put it in data frame
colnames(mutposdf) <- "Mut_Position" 
APCdata <- cbind(APCdata, mutposdf)

### Select the truncating mutations
sgain <- APCdata[APCdata$One_Consequence == "stop_gained", ]
frvar <-APCdata[APCdata$One_Consequence == "frameshift_variant", ]
trunc <- rbind(sgain, frvar)

### KDE plot of APC trucating mutations
ptrunc <- ggplot(trunc, aes(x=Mut_Position)) +  geom_density(aes(y=(..count..)*50), adjust = 1/6)+
  geom_point(stat = "count", colour = "blue")+
  labs(title ="APC KDE for Colorectal cancer", x = "AA position on APC", y = "# of obserced mutations") +
  xlim(0, 2843) +
  theme(plot.title = element_text(size=24))+
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=22), axis.title=element_text(size=22))+ 
  theme(legend.position = c(0.8, 0.8)) +
  geom_vline(xintercept = c(370, 715, 1005, 1235, 1587), color = "red", size=0.5)  + 
  annotate(geom="text", x=80, y=35, label="Zone I", color="black", size = 4) +
  annotate(geom="text", x=500, y=35, label="Zone II", color="black", size = 4) +
  annotate(geom="text", x=800, y=35, label="Zone III", color="black", size = 4) +
  annotate(geom="text", x=1130, y=35, label="Zone IV", color="black", size = 4) +
  annotate(geom="text", x=1400, y=35, label="Zone V", color="black", size = 4) +
  annotate(geom="text", x=2200, y=35, label="Discarded area", color="black", size = 4) 

ptrunc + theme_classic() + # Classic theme
  theme(legend.position = "none") # remove legend


### Sample counts 
SamFreq<- as.data.frame(table(trunc$Tumor_Sample_Barcode))


trunc1 <-  SamFreq %>% filter(Freq == "1")
trunc1data <- trunc[trunc$Tumor_Sample_Barcode %in% trunc1$Var1,]
trunc1data$Mut_Allele <- ifelse(trunc1data$Reference_Allele==trunc1data$Tumor_Seq_Allele1,1,ifelse(dtrunc1data$Reference_Allele==trunc1data$Tumor_Seq_Allele1,2,NA))

###### filter PCA filtered samples and CNA
# we filter again at the end because here we also have to remove samples from samples with no APC truncation
#trunc1data <- trunc1data %>% filter(trunc1data$Tumor_Sample_Barcode != "TCGA-A6-5656-01" & trunc1data$Tumor_Sample_Barcode != "TCGA-A6-2684-01" & trunc1data$Tumor_Sample_Barcode != "TCGA-AZ-4323-01")


#CNV <- read.csv("./produced_data/APC_CNV") 
#colnames(CNV) <- c("cna", "Tumor_Sample_Barcode")
#CNV$Tumor_Sample_Barcode <- sapply(CNV$Tumor_Sample_Barcode, function(x){substring(x, 1, 15)})
#a<- CNV %>% filter(cna != 0)
#trunc1data <- trunc1data[!trunc1data$Tumor_Sample_Barcode %in% a$Tumor_Sample_Barcode,]

###########

### KDE plot of APC trucating mutations
ptrunc1 <- ggplot(trunc1data, aes(x=Mut_Position)) +  geom_density(aes(y=(..count..)*50), adjust = 1/4)+
  geom_point(stat = "count", colour = "blue")+
  labs(title ="Kernel Density Estimates for APC", x = "AA position on APC", y = "# of obserced mutations") +
  xlim(0, 2843) +
  ylim(0, 30)+
  theme(plot.title = element_text(size=24))+
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=22), axis.title=element_text(size=22))+ 
  theme(legend.position = c(0.8, 0.8)) +
  geom_vline(xintercept = c(1235, 1587), color = "red", size=0.5)  + 
#  annotate(geom="text", x=80, y=35, label="Zone I", color="black", size = 4) +
#  annotate(geom="text", x=500, y=35, label="Zone II", color="black", size = 4) +
  annotate(geom="text", x=700, y=30, label="Zone I", color="black", size = 5) +
  annotate(geom="text", x=700, y=28, label="n = 83", color="black", size = 4) +
 # annotate(geom="text", x=1130, y=35, label="Zone IV", color="black", size = 4) +
  annotate(geom="text", x=1412, y=30, label="Zone II", color="black", size = 5) +
  annotate(geom="text", x=1400, y=28, label="n = 103", color="black", size = 4) +
  annotate(geom="text", x=2200, y=30, label="Zone III", color="black", size = 5) +
  annotate(geom="text", x=2200, y=28, label="n = 1", color="black", size = 4) 

ptrunc1 + theme_classic() + # Classic theme
  theme(legend.position = "none") # remove legend

#### ##
### Remove sample in Discarted area from the data
trunc1data<- trunc1data %>% filter(!Tumor_Sample_Barcode == "TCGA-CM-5860-01")

## Part 4. Tables with samples in zones for samples with 

# for 1 truncating mutation
trunc1data$Zones <- sapply(trunc1data$Mut_Position, function(x){if (x<1235) {1} else{2}}) # else if ( 1235 < x && x < 1587) {2} else {"Discarted area"}})
n1tr<- as.data.frame(table(trunc1data$Zones))
n1tr


########### start with sammples with only 1 truncating mutation
##### Place samples in classes 
# Class 0: samples with no mutated APC and samples with mutated but not truncated
temp<- unique(trunc$Tumor_Sample_Barcode) # We excluded sample "TCGA-CM-5860-01"
temp2<- matchedsamples[!matchedsamples %in% temp]
Class_0_samples<- as.data.frame(temp2[!temp2 %in% "TCGA-CM-5860-01"])# samplein discarded area
Class_0_samples$Class <- 0
colnames(Class_0_samples) <- c("Tumor_Sample_Barcode", "Class")

# class 1-5 : samples with only One truncating mutation
trunc1data$Zones
trunc1data$Class <- sapply(trunc1data$Zones, function(x){if (x==1) {1}  else {2}})
trunc1data$Class
Class_1to5_samples <-  trunc1data %>% select(Tumor_Sample_Barcode, Class)

# merge the classes
samp_to_class <- rbind(Class_0_samples, Class_1to5_samples)
samp_to_class%>% group_by(Class) %>% count()
samp_to_class # 385 samples
table(samp_to_class$Class)

# remove PCA filtered samples
#samp_to_class <- samp_to_class %>% filter(samp_to_class$Tumor_Sample_Barcode != "TCGA-A6-5656-01" & samp_to_class$Tumor_Sample_Barcode != "TCGA-A6-2684-01" & samp_to_class$Tumor_Sample_Barcode != "TCGA-AZ-4323-01")

#save samples and classes as a csv
write.csv(samp_to_class,"1_trunc_samples.prefiltered", row.names = FALSE)


####
# look at KRAS in our classes
KRAS <- SNVpost %>% filter(Hugo_Symbol == "KRAS")
a<- unique(KRAS$Tumor_Sample_Barcode)
samp_to_class <- read.csv("./1_trunc_samples.filtered")
c0 <- samp_to_class %>% filter(Class == 0)
duplicated(c0)
c1 <- samp_to_class %>% filter(Class == 1)
c2 <- samp_to_class %>% filter(Class == 2)

sum(c0$Tumor_Sample_Barcode %in% a)
sum(c1$Tumor_Sample_Barcode %in% a)
sum(c2$Tumor_Sample_Barcode %in% a)

chigh <- KRAS %>% filter(IMPACT == "HIGH")
a<- unique(chigh$Tumor_Sample_Barcode)

