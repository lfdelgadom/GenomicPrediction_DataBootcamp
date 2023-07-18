###########################################
# Demo on marker data loading, filtering, imputation and basics of model fitting and testing factors affecting accuracy.                      
# prof: Aaron Lorenz, lore0149@umn.edu 
# Student: Luis Fernando Delgado Munoz, luis.delgado@cgiar.org
# Date: July, 2023
###########################################

# Start with a clean working environment by removing all objects stored in R's memory
rm(list=ls(all=TRUE))

# Source in functions to be used. This also loads in all the libraries you need. Change the path to point R towards where you saved this script.
source("D:\\OneDrive - CGIAR\\Data Analysis\\GenomicPrediction_DataBootcamp\\script\\bootcamp_functions.R")


# Set directory where you have the data files as working dir.. 
workDir <- "D:\\OneDrive - CGIAR\\Data Analysis\\GenomicPrediction_DataBootcamp\\data\\"
setwd(workDir)


# Import phenotypic data
pheno_data <- read.csv("soyname_pheno_all.csv")

# summary pheno data
summary(pheno_data)

# RIL = Recombinant inbred lines

# Load in vcf file containing genotype data. 
# Use package vcfR to read in and work with vcf file.
infileVCF <- "soynam_geno_sub.vcf"
vcf <- read.vcfR(infileVCF, verbose = FALSE)


# Set type of imputation method for use below, either naive or markov

impMethod <- "naive"
# impMethod <- "markov"

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

# Add missing data
# Calculating the total number of elements in a matrix
totEelm <- dim(geno_num5)[1]*dim(geno_num5)[2]

# Select a seed number
set.seed(123)

# sample.int() is a function that generates random samples from a given set of 
# integers. Here, it is used to generate a vector of 10,000 randomly selected 
# indices from the total number of elements in the matrix (totEelm).
addNA <- sample.int(n = totEelm, size = 10000) 


# as.matrix(geno_num5) converts geno_num5 into a matrix object 
# (if it was not already in matrix format).
geno_num5 <- as.matrix(geno_num5)

# the selected indices (addNA) are used to assign NA values to 
# the corresponding elements in the matrix 
geno_num5[addNA] <- NA


# Impute genotype data using either naive imputation or Markov chain implemented in the NAM package
# If impMethod is equal to "naive", then replaces missing values (NA) in the 
# matrix geno_num5 with the mean value of the respective column. 
if (impMethod == "naive") geno_imp <- replaceNAwithMean(geno_num5)

#markov() is a function that performs imputation using a Markov chain-based method. 
if (impMethod == "markov") geno_imp <- markov(apply(geno_num5, 2, as.numeric))

# setting rownames and colnames  
rownames(geno_imp) <- rownames(geno_num5)
colnames(geno_imp) <- colnames(geno_num5)


# Here, you can reduce the number of individuals in the dataset to test the 
# effect of training population size if you wish.

# Select a seed number
set.seed(1234)

# here we can change the training population individual numbers
ssNdx <- sample.int(n = dim(pheno_data3)[1], size=500) 


geno_imp_sub <- geno_imp[ssNdx, ]
pheno_sub <- pheno_data3[ssNdx, ]

rownames(geno_imp_sub) <- pheno_sub$RIL
### Fit some genomic prediction models to the data ----


# Fit an RR-BLUP model using the rrBLUP package for seedsize trait
# Calculates maximum-likelihood (ML/REML) solutions for mixed models
# extracting the estimated marker effects
rrModel <- mixed.solve(y = pheno_sub$Seedsize, Z=geno_imp_sub)
names(rrModel)

# rrModel$u retrieves the estimated random effects from the rrModel object
mrk_effs_RR <- rrModel$u

# Use marker effects to calculate genomic estimated breeding values of 
# individuals in the training set by using. Here we extract the intercept. 
# We add it back.

# rrModel$beta retrieves the estimated fixed effects from the rrModel object, 
#which was obtained by fitting the mixed model earlier.
int <- as.numeric(rrModel$beta)

# %*% is the matrix multiplication operator in R. It performs matrix 
# multiplication between geno_imp_sub and mrk_effs_RR, which represents 
# the estimated marker effects obtained from the mixed model.
# calculating the genomic estimated breeding values (GEBVs) based on a mixed model analysis. 
gebv_rr <- int + geno_imp_sub%*%mrk_effs_RR


# Calculating a genomic relationship matrix using rrBLUP package and fitting a G-BLUP model
# The GRM quantifies the genetic similarity between individuals based on their genomic information
G <- A.mat(geno_imp_sub)
heatmap(G)

# Gblup model for seedsize trait
gblupModel <- kin.blup(data=pheno_sub, geno='RIL', pheno='Seedsize', K=G)
names(gblupModel)

# GBLUP solution for the genetic values
gblupGebv <- gblupModel$g

# Correlation among genomic estimated breeding values and phenotipic values
cor(gebv_rr, pheno_sub$Seedsize)

# Compare GEBVs from ridge regression BLUP to G-BLUP
cor(gebv_rr, gblupGebv)
plot(gebv_rr, gblupGebv)



# Cross-validation analysis -----------------------------------------------
# Now extend this to perform a 10-fold cross-validation analysis
# technique used to assess the performance and evaluate the predictive 
# ability of a genomic prediction model.

# This works if my total sample size is divisible by 10. If not, need to subset so it is.
# The code you provided involves shuffling (randomizing) the order of rows in 
# a matrix and applying the same order to another matrix. 
ndxShuf <- sample(1:dim(geno_imp_sub)[1], dim(geno_imp_sub)[1])

# These lines of code use the shuffled indices to reorder the rows in the 
# pheno_sub and geno_imp_sub matrices. 
pheno_shuf <- pheno_sub[ndxShuf, ]
geno_imp_sub_shuf <- geno_imp_sub[ndxShuf, ]

# The floor function is commonly used to round down a number to the nearest integer.
# creates a sequence of numbers from 1 to the rounded-down division result, 
# assigning it to the variable cnt.
cnt <- 1:floor(length(ndxShuf)/10) 

# vector(length = length(ndxShuf)) creates a logical vector of the specified length, 
# filled with FALSE values by default. The resulting vector is assigned to the variable pred_stor.
pred_stor <- vector(length=length(ndxShuf))

for (i in 1:10){
  pheno_trn <- pheno_shuf
  # This step introduces missing values (NAs) in the training set, 
  # where the number of missing values corresponds to one-tenth of the total 
  # number of shuffled indices.
  pheno_trn$Seedsize[cnt] <- NA
  
  
  # A mixed model (rrModel) is fitted using the training data (pheno_trn and geno_imp_sub_shuf).
  rrModel <- mixed.solve(y = pheno_trn$Seedsize, Z = geno_imp_sub_shuf)
  
  # The marker effects (mrkEffsRR) are extracted from the rrModel object.
  mrkEffsRR <- rrModel$u
  
  # The intercept (int) is obtained from rrModel$beta, converted to numeric, 
  # and added to the marker effects (mrkEffsRR) multiplied by the marker genotypes 
  int <- as.numeric(rrModel$beta)
  
  # This step calculates the genomic estimated breeding 
  # values (gebv_rr) for the individuals in the training set.
  gebv_rr <- int + geno_imp_sub_shuf%*%mrkEffsRR
  
  # The calculated genomic estimated breeding values (gebv_rr) 
  # at the corresponding indices in cnt are stored in the pred_stor vector.
  pred_stor[cnt] <- gebv_rr[cnt]
  
  # The cnt vector is updated by adding the floor value of the length of 
  # ndxShuf divided by 10. This ensures that in the next iteration, 
  # the subsequent one-tenth portion of shuffled indices will be selected.
  cnt <- cnt + floor(length(ndxShuf)/10)
}

# Overall, the for loop performs these steps iteratively, each time using a 
# subset of shuffled data for training, estimating marker effects, 
# calculating genomic estimated breeding values, and storing the predictions 
# in pred_stor. The loop incrementally progresses through the shuffled indices, 
# allowing the analysis to be performed on different training sets and enabling 
# evaluation of the model's performance.

# By computing the correlation, this code provides an assessment of how well 
# the predicted values from the genomic prediction model align with the observed phenotypic values.
cor(pred_stor, pheno_shuf$Seedsize)

x11(); plot(pred_stor, pheno_shuf$Seedsize)
hist(pred_stor)










