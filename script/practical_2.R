###########################################
# Demo on using the BGLR package to fit some different genomic prediction models. 
# Note that this script is the same as for Practical 1 up until Line 174, 
# except in this case we will be using the whole SoyNAM population, not just a subset.                      
# Aaron Lorenz, lore0149@umn.edu, June, 2023
###########################################

# Start with a clean working environment by removing all objects stored in R's memory
rm(list=ls(all=TRUE))


# Source in functions to be used. This also loads in all the libraries you need. Change the path to point R towards where you saved this script.
source("D:\\OneDrive - CGIAR\\Data Analysis\\GenomicPrediction_DataBootcamp\\script\\bootcamp_functions.R")


# Set directory where you have the data files as working dir.. 
workDir <- "D:\\OneDrive - CGIAR\\Data Analysis\\GenomicPrediction_DataBootcamp\\data\\"
setwd(workDir)


# Import phenotypic data: 5528 accessions
pheno_data <- read.csv("soyname_pheno_all.csv")


# Load in vcf file containing genotype data. 
# Use package vcfR to read in and work with vcf file.
# It is necessary download the file from google drive link provided
# https://drive.google.com/file/d/1hMF6mfjvc7VCCXNH9JvjFSmAKYHAdB8M/view?usp=drive_link
infileVCF <- "SoyNAM_Geno.vcf"
vcf <- read.vcfR(infileVCF, verbose = FALSE)


# -------------------------------------------------------------------------

# Set type of imputation method for use below, either naive or markov
impMethod <- "naive"

# Format, manipulate and filter genotype data ----

# Converting VCF file format to numerical matrix format that 
# can be fit in statistical models
geno_num <- convert_vcf_to_genotype(vcf = vcf)

# Match phenotypic data to marker data by merging, 
# then taking back apart again. joined data is a dataframe
joined_data <- left_join(geno_num, pheno_data, by='RIL')

# Select only geno_data. This data frame does not have RIL
geno_num2 <- joined_data[, 1:4292]

# Select only pheno data since 4293 column
pheno_data2 <- joined_data[, 4293:4299]

# -----------------------------------------------------------
# This code defines a custom function called miss.rm. 
# The purpose of this function is to calculate the number 
# of missing values (NA) in a given vector x. It uses the 
# is.na() function to identify NA values, 
# the which() function to determine their positions, 
# and finally, the length() function to count the number 
# of NA values.

miss.rm <- function(x){length(which(is.na(x)))}

# This line of code applies the miss.rm function to a matrix called geno_num2 using the apply() function.
# The MARGIN=2 argument indicates that the function should be applied to each column of the matrix.
# FUN=miss.rm specifies the function to be used.
# The result of the apply() function is then divided by the number of rows in geno_num2 (given by dim(geno_num2)[1]).

# The overall purpose of this code snippet is to calculate the proportion 
# of missing values per column (genomic position) in the geno_num2 matrix. 
# The resulting values will be stored in the mrkNa vector, where each element 
# corresponds to a column in the matrix.

# This filtering action can help identify genomic positions with a high 
# percentage of missing data, which can be useful for subsequent analyses or 
# quality control steps in bioinformatics research.

mrkNa <- (apply(geno_num2, MARGIN = 2, FUN = miss.rm))/dim(geno_num2)[1]


# Finding columns with more than 50% missing values
# identifies the indices of columns (genomic positions) in the 
# mrkNa vector that have a proportion of missing values (NA) greater than 0.5
ndx <- which(mrkNa > 0.5)

# In this case, the columns specified by the indices in ndx are removed using 
# negative indexing (-ndx). The resulting filtered matrix is stored in geno_num3.
if (length(ndx)>0) {
  geno_num3 <- geno_num2[, -ndx]
  
  # If there are no columns with more than 50% missing values, the geno_num2 matrix 
  # remains unchanged, and geno_num3 is assigned the same value as geno_num2.  
} else { 
  geno_num3 <- geno_num2
}


# Filtering out individuals on proportion of missing data
# Calculating the proportion of missing values per row
indNa <- (apply(geno_num3, MARGIN = 1, FUN = miss.rm))/dim(geno_num3)[2]
# Finding rows with more than 50% missing values:
ndx2 <- which(indNa > 0.5)

# If there are indices in ndx2 (i.e., length(ndx2) > 0), 
# it means there are rows with more than 50% missing values
# In this case, the rows specified by the indices in ndx2 are 
# removed using negative indexing (-ndx2).

if (length(ndx2)>0) {
  geno_num4 <- geno_num3[-ndx2, ] 
  
  # If there are no indices in ndx2, meaning there are no rows with more 
  # than 50% missing values, the geno_num3 matrix remains unchanged, 
  # and geno_num4 is assigned the same value as geno_num3.
} else {
  geno_num4 <- geno_num3
}

# If there are indices in ndx2 (i.e., length(ndx2) > 0), 
# it means there are rows (genotypes) with more than 50% missing values
# In this case, the rows specified by the indices in ndx2 are 
# removed from pheno_data using negative indexing (-ndx2)

if (length(ndx2)>0) {
  pheno_data3 <- pheno_data2[-ndx2, ] 
} else {
  pheno_data3 <- pheno_data2
}


# Filter markers based on MAF. Since the genotypes are coded as 1 and 0, 
# the mean of the scores for a given marker is equal to the proportion of 1's. 
# This is a convenient way of doing it, can be modified according to other 
# numerical coding, and factors in imputed marker scores that are not integers.
# Minor Allele Frequency (MAF) refers to the frequency of the less common 
# allele at a given marker in a population.

maf <- apply(geno_num4, MARGIN = 2, FUN = mean, na.rm=T)

# Finding markers with means less than 0.05 & more than 0.95 Minor Allele Frequency 
ndx3 <- which(maf < 0.05 | maf > 0.95) 

if (length(ndx3)>0) {
  geno_num5 <- geno_num4[, -ndx3] 
} else {
  geno_num5 <- geno_num4
}


# Impute genotype data using either naive imputation or Markov chain implemented in the NAM package
# If impMethod is equal to "naive", then replaces missing values (NA) in the 
# matrix geno_num5 with the mean value of the respective column. 
if (impMethod == "naive") geno_imp <- replaceNAwithMean(geno_num5)

#markov() is a function that performs imputation using a Markov chain-based method. 
if (impMethod == "markov") geno_imp <- markov(apply(geno_num5, 2, as.numeric))

# setting rownames and colnames  
rownames(geno_imp) <- rownames(geno_num5)
colnames(geno_imp) <- colnames(geno_num5)


# here we can change the training population individual numbers
# Select a seed number
set.seed(1234)

ssNdx <- sample.int(n = dim(pheno_data3)[1], size=1000) 


geno_imp_sub <- geno_imp[ssNdx, ]
pheno_sub <- pheno_data3[ssNdx, ]

rownames(geno_imp_sub) <- pheno_sub$RIL

############################
# New code for practical 2 #
############################

### BGLR model fitting ----
# Use the BGLR package to fit various types of models. 
# BRR = Bayesian ridge regression, 
# BL = Bayes LASSO, BayesA, BayesB, BayesC 

# Remove some data to perform a validation analysis
# Use line coding to identify RILs by family 

# "...$": This is a regular expression pattern that matches the last three 
# characters of a string (... represents any three characters).
# "": This is the replacement string, which means that the matched pattern 
# (last three characters) will be replaced with nothing (in essence, they will be removed).
fam <- gsub("...$", "", rownames(geno_imp_sub))
ndxFam <- which(fam=="DS11-64")

pheno_sub_trn <- pheno_sub

pheno_sub_trn$Seedsize[ndxFam] <- NA # removing family DS11-64  

G <- A.mat(geno_imp_sub) # Calculate a genomic relationship matrix using the rrBLUP function A.mat

# Save the relationship matrix plot into a pdf
pdf(file = paste0(here::here(), "/", "G_matrix.pdf"))
    
heatmap(G, scale="column", col = terrain.colors(256),
                  main = "Genomic Relationship Matrix Heatmap")
dev.off()



# G matrix scaled
gscale <- scale(G, center = T, scale = T)

# sometimes is usefful work with G matrix scaled and centered.

D <- (as.matrix(dist(gscale, method = "euclidean"))^2)/dim(geno_imp)[2]

h <- 0.5

K <- exp(-h*D)

# The BGLR package has us first create a list specifying the parameters and 
# some inputs such as marker data and model type
ETA <- list(list(K=NULL, X=geno_imp_sub, model='BayesC', probIn=.10)) 

# We can add Genomic relationship matrix scaled as a imput to the model
# When we add the G matrix, it is mandatory remove the genotypic data, 
# because G Matrix already contents Geno infomation.

# Fitting a Single Kernel Model in BGLR
# The BGLR package has us first create a list specifying the parameters and 
# some inputs such as marker data and model type
ETA <- list(list(K=K, X=NULL, model='RKHS', probIn=.10)) 

model_bglr <- BGLR(y=pheno_sub_trn$Seedsize, ETA=ETA, burnIn=500, nIter=2000, 
                   saveAt='RKHS_h=0.5_') 

# We can add Genomic relationship matrix as a imput to the model

"ETA <- list(list(K=G, X=geno_imp_sub, model='BayesC', probIn=.10)) 
#The BGLR package has us first create a list specifying the parameters and 
some inputs such as marker data and model type"


#model_bglr <- BGLR(y=pheno_sub_trn$Seedsize, ETA=ETA, burnIn=500, nIter=2000, verbose=TRUE) 

# BGLR directly outputs the genotype predictions as yHat
gebv_bglr <- model_bglr$yHat  

# Extract marker effect predictions from model object. 
# Try different models, changing the name of the object storing the effect 
# (e.g,, "bhat_brr") and plot them against one another on a scatter plot.
bhat <- model_bglr$ETA[[1]]$b

# Here is a way to make a trace plot
plot(bhat^2, ylab='Estimated squared marker effect', type='o')


##Correlate predictions of RILs left out of the analysis, with predictions
cor(pheno_sub$Seedsize[ndxFam],  gebv_bglr[ndxFam])

family_DS11_64 <- 
  tibble(pheno_sub = c(pheno_sub$Seedsize[ndxFam]),
       gebv_bglr = c(gebv_bglr[ndxFam]))

library(ggplot2)
library(ggpubr)

  family_DS11_64 %>% ggplot(aes(x = pheno_sub, y = gebv_bglr)) +
  geom_point(shape = 19, 
             size = 3) + 
    stat_regline_equation() + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    theme_bw()
  
ggsave("RKHS.png", units = "in", dpi = 300, width = 6, height = 6)
  
# -------------------------------------------------------------------------

# Fit a multi-kernel model using BGLR to treat some large-effect QTL as fixed effects, 
# and remaining QTL as random effects. QTL here were previously declared 
# significant using a GWAS analysis. SNP positions of QTL were 1926, 829, 683, 678.

qtl <- c(1926, 829, 683, 678)

ETA_mk <- list(list(X=geno_imp_sub[, qtl], model='FIXED', probIn=.10), 
               list(K=G, X=geno_imp_sub[, -qtl], model='RKHS', probIn=.10))

model_bglr_mk <- BGLR(y=pheno_sub_trn$Seedsize, ETA=ETA_mk, burnIn=500, nIter=2000, 
                      verbose=FALSE)

gebv_bglr_mk <- model_bglr_mk$yHat

cor(pheno_sub$Seedsize[ndxFam],  gebv_bglr_mk[ndxFam])
plot(pheno_sub$Seedsize[ndxFam],  gebv_bglr_mk[ndxFam])

