# rm(list = ls())
# gc()

args = commandArgs(trailingOnly=TRUE)
print(args)


dir_path = args[1]
setwd(dir_path)
source("../src.R")


##Load Data
data = read.csv(args[2], row.names = 1)
metadata = read.csv(args[3], row.names = 1)

#set parameters
multi_GRM_flag = 0 #1 = multiple variance components; 0 = otherwise 
fixed_effect_flag = 0 #1 = the model incorporates fixed effects; 0 = otherwise 

#Impute metadata
EM_flag = 1
if(EM_flag == 1){
  
  library("Amelia")
  
  a.out <- amelia(metadata, m = 100, cs = "Depth_m" ,p2s = 2)
  I_metadata = a.out$imputations$imp100
  
  I_metadata_New = na.omit(I_metadata)
  
  for(j in 1:dim(I_metadata_New)[2]){
    
    min_tmp = min(I_metadata_New[,j])
    Num_Neg_values = which(I_metadata_New[,j] < 0)
    if(length(Num_Neg_values) > 0)
      I_metadata_New[,j] = I_metadata_New[,j] + abs(min_tmp)
  }
  
  metadata = I_metadata_New
  
}


# Extract only those samples in common between the two tables - 712 samples
common.sample.ids <- intersect(rownames(metadata), rownames(data))
data <- data[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# Double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

index <- c()
for(j in 1:dim(data)[2]){
  
  split = strsplit(colnames(data)[j], split = "eggNOG_OGs_")
  tmp = split[[1]][2]
  index = c(index, tmp)
}


X = as.matrix(data)
n = length(which(is.na(X)))

#Imputing missing values using uniform distribution in [0,0.1]
X[which(is.na(X))] = runif(n, min = 0, max = 10^-6)

##Create the kinship matrix - using all the eggNOG genes (after imputation)

dat = metadata
scaled.dat <- scale(dat) #normalize each column
scaled.dat[is.na(scaled.dat)] = 0

#Create the kinship matrix
Kinship_test = cosine(scaled.dat)
ev <- eigen(Kinship_test)
eval <- ev$values

#in case the kinship matrix is not positive semidefinite
if(any(eval) < 0){
  
  nc.  <- nearPD(Kinship_test, conv.tol = 1e-7) # default
  Kinship_test <- nc.$mat 
}


###STEP 3 - Calculating 'Time_explainability'###

#one variance component
if(multi_GRM_flag == 1){
  config_flag = "2_GRM"
  TE_path = "eggNOG_explainability"
}else{
  config_flag = "1_GRM"
  TE_path = "eggNOG_explainability"
}


if(dir.exists(file.path(paste0(dir_path, "/", TE_path)))){
  
  VE_Results = Calculate_VE_func(path_hsq_files = paste0(dir_path,"/", TE_path),
                                 Data = scaled.dat, All_individuals = dim(scaled.dat)[1],
                                 num_time_points = 1,Fixed_effect_flag = fixed_effect_flag ,ind_1 = NULL,
                                 config = config_flag, hsq_txt_tmp = "GRM_reml__")
  
  colnames(VE_Results) = c("variance_explained_Env", "SD_variance_explained_Env", 
                           "p_value", 
                           "logL0", "logL", "gene_index")
  
  
  VE_Results = VE_Results[order(VE_Results$gene_index),]
  p_adjust = p.adjust(p = VE_Results$p_value, method = "BH", n = length(VE_Results$p_value))
  
  VE_Results$p_value_adjusted = p_adjust
  
  VE_Results = VE_Results[,c("variance_explained_Env", "SD_variance_explained_Env",
                             "p_value", 
                             "logL0", "logL", "gene_index")]
  

  VE_Results$eggNOG = index[VE_Results$gene_index]
  

  print(VE_Results)
  
}

VE_Results_sort <- VE_Results[order(-VE_Results$variance_explained_Env),]
par(mfrow = c(1,1))
boxplot(head(VE_Results_sort$rank, 100), tail(VE_Results_sort$rank, 100), names = c("Env Variance explained Top",
                                                                                    "Env Variance explained Bottom"),
        ylab = "Gene expression rank", main = "p-value < 5e-07 (Wilcoxon & t-test)")

wilcox.test(head(VE_Results_sort$rank, 100),
            tail(VE_Results_sort$rank, 100))

t.test(head(VE_Results_sort$rank, 100),
            tail(VE_Results_sort$rank, 100))

setwd(dir_path)
dir.create("Results")
write.csv(VE_Results, file = "VE_Results_eggNOG_Env_kinship.csv")

