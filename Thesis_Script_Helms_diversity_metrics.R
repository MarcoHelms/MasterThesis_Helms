# Calculation of taxonomic and functional diversity metrics on an alpha and beta
# diversity level

# 1. Create and match data sets ---------------------------------------------------
# Upload a trait value data sets in the following dimensions
# Rows      = Taxa (names)
# Columns   = Functional traits (names)
# Fill      = Values per taxa per trait
spiders_trt <- read.table("spiders_coded_traits.txt", header=T, row.names = 1, sep = "\t")


# Upload a taxa abundance data sets in the following dimensions
# Rows      = Sample sites (names)
# Columns   = Taxa (names)
# Fill      = Abundance per taxa per site
spiders_tax <- read.table("Garden_spiders_matrix.txt", header=T, row.names = 1, sep = "\t")

# Substitute symbols in data that are not matching
rownames(spiders_trt) <- gsub('_', '.',rownames(spiders_trt))

# Match the data by sampled taxa names, order them and fit them together
spiders_trt <- spiders_trt[rownames(spiders_trt) %in% colnames(spiders_tax), ]
spiders_tax <- spiders_tax[,colnames(spiders_tax) %in% rownames(spiders_trt), ]

spiders_tax <- spiders_tax[,order(colnames(spiders_tax))]
cbind(colnames(spiders_tax), rownames(spiders_trt))


# Check if the sampled taxa are identical in both data sets share dimansions
identical(colnames(spiders_tax), rownames(spiders_trt))
dim(spiders_tax); dim(spiders_trt)


# 2. Calculate Alpha Taxonomic Metrics --------------------------------------------
require(vegan)

# Prepare data for analyses
trt = spiders_trt
sp = spiders_tax
rownames(sp) <- paste("site", 1:nrow(sp), sep="_")

# Calculate richness and diversity
rich <- specnumber(sp)
sp = sp[rich>=3,]

rich = specnumber(sp)
sha = vegan::diversity(sp)

# 3. Calculate Community weighted means for functional traits----------------------
CWM_gen <- function(L,T, Chessel = TRUE){
  
  T<-as.matrix(T)
  L<-as.matrix(L)}
standardize_w <- function(X,w){
  ones <- rep(1,length(w))
  Xc <- X - ones %*% t(w)%*% X
  Xc / ones%*%sqrt(t(ones)%*%(Xc*Xc*w))} 

CWM_calc <- function(trt, site){
  CWM = matrix(NA, nrow(site), ncol(trt))
  colnames(CWM) = colnames(trt)
  rownames(CWM) = rownames(site)
  for(i in 1:ncol(trt)) CWM[rowSums(sp[,!(is.na(trt[,i]))])>0,][,i] <- CWM_gen(L=sp[,!(is.na(trt[,i]))],T=na.omit(trt[,i]))[,1]
  return(CWM)}


CWM = CWM_calc(spiders_trt, site=sp)

for(i in 1:ncol(CWM)){
  print(i)
  CWM[,i][is.na(CWM[,i])]=mean(na.omit(CWM[,i]))
}

# Create a new data frame with the Community weighted means
spiders_traits_CWM <- data.frame(CWM)


# 4. Calculate Alpha Functional Metrics--------------------------------------------
library(caret)
library(recipes)
require(caret)
require(vegan)

# Prepare data and reduce dimesionality
trt.mis.model = preProcess(trt, "knnImpute")
trt.mis.pred = predict(trt.mis.model, trt); head(trt.mis.model)

trt.pca <- prcomp(trt.mis.pred, scale. = T, center = T)
cumsum(trt.pca$sdev/sum(trt.pca$sdev))
trt.scaled <- trt.pca$x[,1:2]

# Calculate Functional Dispersion (Laliberté 2010)
calc_mFD <- function(trt, sp){
  
  require("mFD")
  
  trt.cat = data.frame(trait_name=colnames(trt), 
                       trait_type = rep("Q", ncol(trt)),                       
                       trait_type = rep(1, ncol(trt)),
                       fuzzy_name =NA)
  
  
  sp_dist_trt <- mFD::funct.dist(
    sp_tr         = trt,
    tr_cat        = trt.cat,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  
  FD2max <- mFD::alpha.fd.hill(
    asb_sp_w = as.matrix(sp), 
    sp_dist  = sp_dist_trt, 
    tau      = "max", 
    q        = 2)
  
  fspaces <- mFD::quality.fspaces(
    sp_dist             = sp_dist_trt,
    maxdim_pcoa         = 10,
    deviation_weighting = "absolute",
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  sel.axis = which(fspaces$"quality_fspaces" == min(fspaces$"quality_fspaces"))
  sel.axis=2
  fspace.sel <- fspaces$details_fspaces$sp_pc_coord[,1:sel.axis]
  
  
  alpha_fd <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = fspace.sel,
    asb_sp_w         = as.matrix(sp[rich>sel.axis,]),
    ind_vect         = c("fdis", "feve", "fric"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  dim(alpha_fd$functional_diversity_indices)
  dim(sp)
  
  FD_mat_Magneville <- as.data.frame(matrix(NA, 
                                            ncol=ncol(alpha_fd$functional_diversity_indices), 
                                            nrow=nrow(sp)
  ))
  colnames(FD_mat_Magneville) = colnames(alpha_fd$functional_diversity_indices)
  rownames(FD_mat_Magneville) = rownames(sp)
  
  for(i in 1:nrow(sp)){
    if(rownames(FD_mat_Magneville)[i] %in% rownames(alpha_fd$functional_diversity_indices)){
      row.sel = which(rownames(alpha_fd$functional_diversity_indices) == rownames(FD_mat_Magneville)[i])
      FD_mat_Magneville[i,] = alpha_fd$functional_diversity_indices[row.sel,]  
    }else FD_mat_Magneville[i,] = rep(NA, ncol(FD_mat_Magneville))
  }
  
  return(data.frame(FD_mat_Magneville, Hill_2_max = FD2max$asb_FD_Hill))
}


FD_mFD <- calc_mFD(trt=trt.mis.pred, sp=sp)

# Calculate functional Evenness (Ricotta 2014)
fundist_spe <- function(trt, sp){
  require("mFD")
  
  trt.cat = data.frame(trait_name=colnames(trt), 
                       trait_type = rep("Q", ncol(trt)),                       
                       trait_type = rep(1, ncol(trt)),
                       fuzzy_name =NA
  )
  
  
  sp_dist_trt <- mFD::funct.dist(
    sp_tr         = trt,
    tr_cat        = trt.cat,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  return(sp_dist_trt)  
}

trt.dist <- fundist_spe(trt=trt.mis.pred, sp=sp)
dis <- as.matrix(trt.dist)


FeveR<-function(abundances, distances)
{
  rel_abundance<-sweep(abundances, 1, rowSums(abundances), "/")
  rel_abu_matrix<-as.matrix(rel_abundance)
  n_plot<-nrow(rel_abu_matrix)
  n_species<-ncol(rel_abu_matrix)
  dista_matrix<-as.matrix(distances)
  index_array<-rep(NA, n_plot)
  for(i in 1:n_plot)
  {mat_develop<-rep(rel_abu_matrix[i,],n_species)
  matrix_two<-(t(matrix(mat_develop, nrow=n_species,ncol=n_species)))
  matrix_non<-matrix(mat_develop, nrow=n_species,ncol=n_species)
  uni<-matrix(1,nrow=n_species,ncol=n_species)
  subtraction<- uni-matrix_non
  division<-matrix_two/subtraction
  moltiplication<-division*dista_matrix
  row_sums<-rowSums(moltiplication)
  per_abundance<-row_sums*rel_abu_matrix[i,]
  great_sum<-sum(per_abundance)
  divised<-per_abundance/great_sum
  espress<-which(per_abundance>0)
  S<-length(espress)
  bulla<-rep(NA, length(divised))
  for(l in 1:length(divised))
  {
    bulla[l]<-min(divised[l], 1/S)
  }
  indice<-sum(bulla)
  norm_index<-(indice-(1/S))/(1-(1/S))
  index_array[i]<-norm_index
  }
  return(index_array)
}

FEve.Ricotta <- FeveR(sp, trt.dist)

# Create a data frame of the results
df.FD <- data.frame(richness=rich,shannon = sha,
                    CWM, 
                    FD_mFD$fdis, FEve.Ricotta)

# 5. Calculate Taxonomic Beta Diversity (LCBD)-------------------------------------
library(betapart)
require(adespatial)

LCBD.tax.mult <- beta.multi.abund(sp, index.family = "bray")

out3 = beta.pair.abund(sp, index.family = "bray")
LCBD.tax <- data.frame(turnover=LCBD.comp(out3$beta.bray.bal, sqrt.D=TRUE)[[2]], 
                       nested = LCBD.comp(out3$beta.bray.gra, sqrt.D=TRUE)[[2]], 
                       all = LCBD.comp(out3$beta.bray, sqrt.D=TRUE)[[2]])

# 6. Calculate Functional Beta Diversity (LCBD)-------------------------------------
test.pair<-functional.beta.pair(x=decostand(sp, "pa"), traits=trt.scaled, index.family = "jaccard")
LCBD.fun <- data.frame(turnover=LCBD.comp(test.pair$funct.beta.jtu, sqrt.D=TRUE)[[2]], 
                       nested = LCBD.comp(test.pair$funct.beta.jne, sqrt.D=TRUE)[[2]], 
                       all = LCBD.comp(test.pair$funct.beta.jac, sqrt.D=TRUE)[[2]])

# Create a data frame of all beta diversity metrics
df.LCBD <- data.frame(LCBD.tax.turn = LCBD.tax$turnover,
                      LCBD.tax.nest = LCBD.tax$nested,
                      LCBD.tax.all = LCBD.tax$all,
                      LCBD.fun.turn = LCBD.fun$turnover,
                      LCBD.fun.nest = LCBD.fun$nested,
                      LCBD.fun.all = LCBD.fun$all)


# 7. Create the final diversity metrics data frame for this taxonomic -------------
Com.mat <- data.frame(Richness = df.FD$richness,
                      Shannon = df.FD$shannon,
                      Functional_dispersion = df.FD$FD_mFD.fdis,
                      Functional_evenness = df.FD$FEve.Ricotta,
                      df.LCBD)


#Save the data frame 
write.table(Com.mat, "spiders_community_metrics.txt", sep="\t")