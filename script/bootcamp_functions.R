##Functions for GS bootcamp course

rm(list = ls())
library(pacman)
pacman::p_load(vcfR,
               tibble,
               rrBLUP,
               BGLR,
               PopVar,
               qtl,
               STPGA,
               NAM,
               sommer,
               MASS,
               SoyNAM,
               dplyr)


#Easy impute, replacing missing with marker mean
replaceNAwithMean <- function(mat){
  replaceNAwithMeanVec <- function(vec){
    mu <- mean(vec, na.rm=TRUE)
    vec[is.na(vec)] <- mu
    return(vec)
  }
  return(apply(mat, 2, replaceNAwithMeanVec))
}


#Function to simulate phenotypic values
simulatePhenotypicValues_Corr <- function(nextGenGenoTable,NTraits,MaxPotential,H,GenoCorr,noLoci,simSeed){
  
  set.seed(125+simSeed)
  genoTable <- nextGenGenoTable
  nIndividuals <-nrow(genoTable)
  nMarkers <- ncol(genoTable) 
  nLoci <- noLoci
  maxPotential <- MaxPotential
  nTraits <- NTraits
  
  
  h1 <- sqrt(H[1])
  h2 <- sqrt(H[2])

  S <- GenoCorr
  G <- mvrnorm(nLoci,mu=rep(0,nTraits),Sigma=(S))
  
  max_a <- maxPotential/noLoci 
  min_a <- -max_a
  
  G_Eff <- G
  
  scaleFactor <- c()
  for(nT in 1:nTraits){
    scaleFactor[nT] <-  (max_a[nT]-min_a[nT])/(max(G[,nT])-min(G[,nT]))
    #G_Eff[,nT] <- ((max_a[nT]-min_a[nT])/2) + (G[,nT]*(scaleFactor[nT]))
  }
  
  G_Eff <- G %*%  diag(scaleFactor)
  
  ##### Extract QTL table ######
  QTL_indices_List <- list() 
  QTL_table_List <- list() 
  
  for(nT in 1:nTraits){
    QTL_indices_List[[nT]] <- sample(nMarkers,noLoci)
    QTL_table_List[[nT]] <- genoTable[,c(QTL_indices_List[[nT]])]
    Geno <- cbind(QTL_table_List[[nT]] %*% G_Eff)
  }
  
  Geno_M <- apply(Geno,2,function(x) x+ (max(x)-min(x)/2))
  #Geno_M <- Geno
  
  #### Pheno
  
  
  nTraits <- NTraits
  varG <- apply(Geno_M,2,var)
  varE <- c()
  for(nT in 1:nTraits){
    varE[nT] <- ((1-H[nT])*varG[nT])/H[nT]
  }
  
  S_E <- diag(varE)
  Pheno <- Geno_M + mvrnorm(nrow(Geno_M),mu=rep(0,nTraits),Sigma=S_E)
  
  
  return(Pheno)
}


# The function convert_vcf_to_genotype takes a vcf object as input and performs 
# the sequence of operations to convert it to a genotype matrix. 
convert_vcf_to_genotype <- function(vcf) {
  
  #  extracts the genotype information from the vcf object using 
  # the extract.gt function.
  gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  
  # converts the vcf object to a tibble using the as_tibble function.
  fix_T <- as_tibble(getFIX(vcf))
  
  # creates an empty matrix gt2 with the same dimensions as gt.
  gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
  
  # assigns the column names from gt to gt2
  colnames(gt2) <- colnames(gt)
  
  # applies the gsub function to replace "1/1" with "1" in each column of gt.
  gt2a <- apply(gt, 2, function(x) gsub("1/1", "1", x))
  
  # uses gsub to replace "0/0" with "0" in gt2a.
  gt2b <- gsub("0[/|]0", "0", gt2a)
  
  # replaces "1/0" or "0/1" with "0.5" in gt2b using gsub.
  gt2c <- gsub("[10][/|][10]", "0.5", gt2b)
  
  # replaces missing values ".", represented as "./.", with "NA" in gt2c.
  gt2d <- gsub("\\.[/|]\\.", "NA", gt2c)
  
  # converts the elements of gt2d to numeric using as.numeric.
  gt2d_num <- apply(gt2d, 2, as.numeric)
  
  # adds row names back to gt2d_num using the original row names from gt2d.
  rownames(gt2d_num) <- rownames(gt2d)
  
  # transposes gt2d_num and converts it to a data frame geno_num.
  geno_num <- as.data.frame(t(gt2d_num))
  
  # adds a column named "RIL" to geno_num containing the row names.
  geno_num$RIL <- rownames(geno_num)
  
  # the function returns the resulting geno_num data frame.
  return(geno_num)
}

