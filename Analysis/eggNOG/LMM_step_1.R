# rm(list = ls())
# gc()
#
# args = commandArgs(trailingOnly=TRUE)
# print(args)


dir_path = "~/Dropbox/2019-OceanMicrobes/LMM/Analysis/Results_eggNOG_Env_kinship/"
setwd(dir_path)
source("~/Dropbox/2019-OceanMicrobes/LMM/Analysis/src.R")


##Load Data

##Gene relative abundance
library(data.table)
library(R.utils)

setwd("~/Dropbox/2019-OceanMicrobes/LMM/Data/Gene_data/pnps/")
data = fread("eggNOG_OGs_4_60_5_50_5.csv.gz") #744 samples x 16426 genes
data = data.frame(data)
rownames(data) = data[,1]
data = data[,-1]

setwd("~/Dropbox/2019-OceanMicrobes/LMM/Data/Metadata/")
metadata = read.csv("oceans_metadata.csv", row.names = 1)

#Label shuffling
metadata_shuffle <- metadata
idx_shuffle <- sample(c(1:dim(metadata)[1]), size = dim(metadata)[1])
rownames(metadata_shuffle) = rownames(metadata)[idx_shuffle]

metadata <- metadata_shuffle

#Impute metadata
EM_flag = 1
if(EM_flag == 1){

  library("Amelia")

  a.out <- amelia(metadata[,-8], m = 100, cs = "Depth_m" ,p2s = 2)
  I_metadata = a.out$imputations$imp100

  # length(which(is.na(I_metadata))) #181 missing values
  # length(which(is.na(metadata))) #794

  I_metadata_New = na.omit(I_metadata)

  for(j in 1:dim(I_metadata_New)[2]){

    min_tmp = min(I_metadata_New[,j])
    Num_Neg_values = which(I_metadata_New[,j] < 0)
    if(length(Num_Neg_values) > 0)
      I_metadata_New[,j] = I_metadata_New[,j] + abs(min_tmp)
  }

  metadata = I_metadata_New

}


# Extract only those samples in common between the two tables - 687 samples
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

num_missing_samples = c()

for(j in 1:dim(data)[2]){

  num_missing_samples[j] = length(which(is.na(data[,j])))
}

setwd(dir_path)
Info_eggNOG_OGs_4_60_5_50_5 = data.frame(colnames(data), num_missing_samples)
write.csv(Info_eggNOG_OGs_4_60_5_50_5, file = "Info_eggNOG_OGs_4_60_5_50_5_all.csv", row.names = F)


# ind_to_remove = which(is.na(data[,k]))
# X = as.matrix(data[-ind_to_remove,-k])
# y = data[,k][-ind_to_remove]

X = as.matrix(data)
n = length(which(is.na(X)))

#Imputing missing values using uniform distribution in [0,0.1]
X[which(is.na(X))] = runif(n, min = 0, max = 0.1)

##Create the kinship matrix - using all the Kegg genes (after imputation)

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

setwd(dir_path)
dir.create("Kinship_mat_files")
setwd("./Kinship_mat_files")

#GRM file for the metadata matrix W_1
GRM_1 = create_GRM(data_OTUs = scaled.dat,
                   Kinship_mat = Kinship_test ,
                   start_t = 1,
                   end_t = dim(Kinship_test)[1],
                   id =  c(1:dim(Kinship_test)[1]),
                   save_flag = 1,
                   save_path = paste(dir_path, "/Kinship_mat_files", sep = ""),
                   GRM_title = "Kinship_env")


dat_2 = X
scaled.dat_2 <- scale(dat_2) #normalize each column
scaled.dat_2[is.na(scaled.dat_2)] = 0

#Create the kinship matrix
Kinship_test_2 = cosine(scaled.dat_2)
ev <- eigen(Kinship_test_2)
eval <- ev$values

#in case the kinship matrix is not positive semidefinite
if(any(eval) < 0){

  nc.  <- nearPD(Kinship_test_2, conv.tol = 1e-7) # default
  Kinship_test_2 <- nc.$mat
}

#GRM file for the gene matrix W2
GRM_2 = create_GRM(data_OTUs = scaled.dat_2,
                   Kinship_mat = Kinship_test_2 ,
                   start_t = 1,
                   end_t = dim(Kinship_test_2)[1],
                   id =  c(1:dim(Kinship_test_2)[1]),
                   save_flag = 1,
                   save_path = paste(dir_path, "/Kinship_mat_files", sep = ""),
                   GRM_title = "Kinship_genes")


create_mgrm_file(t_min = 1, t_max = 2, num_GRM = 1, prediction_flag = 0,
                 saving_path = dir_path)


###Craete Fixed effects files

# if(dir.exists(file.path(paste(dir_path, "/Fixed_effects", sep = ""))) == FALSE){
#   create_Fixed_effect_files(Data = metadata,
#                                     path_to_save = dir_path) }

###Craete 'phenotype' files
if(dir.exists(file.path(paste(dir_path, "/Phen_files", sep = ""))) == FALSE){
  create_pheno_files(Data = scaled.dat_2, ind = c(1:dim(X)[2]), path_to_save = dir_path)}

setwd(dir_path)
write.table(c(1:dim(X)[2]), file = "ind_NEW.txt", row.names = F, col.names = F, quote = F)
