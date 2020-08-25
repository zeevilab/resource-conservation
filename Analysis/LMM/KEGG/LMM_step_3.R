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
  
  split = strsplit(colnames(data)[j], split = "KEGG_ko_ko.")
  tmp = split[[1]][2]
  index = c(index, tmp)
}


##Ranking using gene expression
setwd(dir_path)
RNA = read.csv(args[4], row.names = 1)

#Data collection:
#Seawater was size fractionated in situ onto large fraction (5µm) and small fraction (0.22µm) filters

read_counts_all <- RNA[,c(3:18)]
read_counts_large <- RNA[RNA$fraction == "large",c(3:18)]
read_counts_small <- RNA[RNA$fraction == "small",c(3:18)]

seq_depth = c()
seq_depth_large = c()
seq_depth_small = c()

for(j in 1:dim(read_counts_large)[2]){
  
  seq_depth_large[j] = sum(na.omit(read_counts_large[,j]))
  seq_depth_small[j] = sum(na.omit(read_counts_small[,j]))
  
}

for(j in 1:dim(read_counts_all)[2]){
  
  seq_depth[j] = sum(na.omit(read_counts_all[,j]))
  
}



KO <- as.character(unique(RNA$KO))
KO = KO[KO != ""] #4029

#exp per fraction
mean_exp_large = c()
mean_exp_small = c()

#both fractions combined
mean_exp = c()
max_exp = c()
var_exp = c()

fraction_flag = 1 #0 = both fractions combined; 1 = per fraction

for(i in 1:length(KO)){
  
  
  if(fraction_flag == 1){
    
    tmp = RNA[RNA$KO == KO[i],]
    tmp_large = tmp[tmp$fraction == "large",c(3:18)]
    tmp_small = tmp[tmp$fraction == "small",c(3:18)]
    
    rel_abn_large = matrix(NA, ncol = dim(tmp_large)[2], nrow = dim(tmp_large)[1])
    rel_abn_small = matrix(NA, ncol = dim(tmp_small)[2], nrow = dim(tmp_small)[1])
    
    for(j in 1:dim(tmp_large)[2]){
      
      rel_abn_large[,j] = tmp_large[,j]/seq_depth_large[j]
      rel_abn_small[,j] = tmp_small[,j]/seq_depth_small[j]
      
    }
    
    rel_abn_large[is.na(rel_abn_large)] = 0
    rel_abn_small[is.na(rel_abn_small)] = 0
    
    mean_exp_large[i] = mean(apply(rel_abn_large, 2, sum))
    mean_exp_small[i] = mean(apply(rel_abn_small, 2, sum))
    
  }
  
  if(fraction_flag == 0){
    
    tmp = RNA[RNA$KO == KO[i],]
    tmp_all_fractions = tmp[,c(3:18)]
    
    rel_abn_all_fractions = matrix(NA, ncol = dim(tmp_all_fractions)[2], nrow = dim(tmp_all_fractions)[1])
    
    for(j in 1:dim(tmp_all_fractions)[2]){
      
      rel_abn_all_fractions[,j] = tmp_all_fractions[,j]/seq_depth[j]
      
    }
    
    rel_abn_all_fractions[is.na(rel_abn_all_fractions)] = 0
    
    mean_exp[i] = mean(apply(rel_abn_all_fractions, 2, sum))
    max_exp[i] = quantile(apply(rel_abn_all_fractions, 2, sum), probs = 1) #Q3
    var_exp[i] = var(apply(rel_abn_all_fractions, 2, sum)) #Q3
    
  }
  
}

names(mean_exp_large) = names(mean_exp_small) = KO
KO_intersection = intersect(KO, index) #1764 KOs

# which(index == "K01915")

small_fraction_mean_exp = sort(mean_exp_small[names(mean_exp_small) %in% KO_intersection], decreasing = T)

View(data.frame(names(small_fraction_mean_exp), small_fraction_mean_exp, rank(small_fraction_mean_exp)))

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


###STEP 3 - Calculating 'Time_explainability'###

#one variance component
if(multi_GRM_flag == 1){
  config_flag = "1_GRM_FF"
  TE_path = "Kegg_explainability_FF"
}else{
  config_flag = "1_GRM"
  TE_path = "Kegg_explainability"
}


if(dir.exists(file.path(paste0(dir_path, "/", TE_path)))){
  
  VE_Results = Calculate_VE_func(path_hsq_files = paste0(dir_path,"/", TE_path),
                                 Data = scaled.dat, All_individuals = dim(scaled.dat)[1],
                                 num_time_points = 1,Fixed_effect_flag = fixed_effect_flag ,ind_1 = NULL,
                                 config = config_flag, hsq_txt_tmp = "GRM_reml__")
  
  colnames(VE_Results) = c("variance_explained", "SD_variance_explained", 
                           "p_value", "logL0", "logL", "gene_index")
  
  
  VE_Results = VE_Results[order(VE_Results$gene_index),]
  p_adjust = p.adjust(p = VE_Results$p_value, method = "BH", n = length(VE_Results$p_value))
  
  VE_Results$p_value_adjusted = p_adjust
  
  VE_Results = VE_Results[,c("variance_explained", "SD_variance_explained",  
                             "logL0", "logL","gene_index", "p_value_adjusted")]
  
  VE_Results$KO = index[VE_Results$gene_index]
  
  rank = c()
  for(j in 1:length(VE_Results$KO)){
    
    if(length(which(names(small_fraction_mean_exp) == VE_Results$KO[j])) > 0)
      rank[j] = length(small_fraction_mean_exp) - which(names(small_fraction_mean_exp) == VE_Results$KO[j])+ 1
    else
      rank[j] = NA
  }

  VE_Results$rank = rank
  

  View(VE_Results)
  
}

VE_Results_sort <- VE_Results[order(-VE_Results$variance_explained),]
boxplot(head(VE_Results_sort$rank, 100), tail(VE_Results_sort$rank, 100), names = c("Variance explained Top",
                                                                                    "Variance explained Bottom"),
        ylab = "Gene expression rank", main = "p-value < 1e-09 (Wilcoxon & t-test)")

wilcox.test(head(VE_Results_sort$rank, 100),
            tail(VE_Results_sort$rank, 100))

t.test(head(VE_Results_sort$rank, 100),
            tail(VE_Results_sort$rank, 100))

setwd(dir_path)
dir.create("Results")
write.csv(VE_Results, file = "VE_Results_KEGG_Env_kinship.csv")

