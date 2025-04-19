##### O2
# srun --pty -p interactive -t 0-12:00 --mem 40G /bin/bash

library(dplyr)
library(data.table)
library(kit)

###########
# FUNCTIONS
retrieve.cos.sim <- function(pair) {
  cos_sim_matrix[[pair[1], pair[2]]]
}

retrieve.cos.sim.2 <- function(matrix, pair) {
  matrix[[pair[1], pair[2]]]
}

# input: U is the embedding matrix (each row is a variable)
cos.sim.matrix <- function(U)
{
  U.norm = U / apply(U, 1, norm,'2')
  return(U.norm %*% t(U.norm))
}
cos.sim.two.matrices <- function(x, y)
{
  x.norm = x / apply(x, 1, norm, '2')
  y.norm = y / apply(y, 1, norm, '2')
  return(x.norm %*% t(y.norm))
}
rescale.max <- function(U)
{
  maxes = apply(U, 1, function(x) max(abs(max(x)), abs(min(x))))
  maxes[maxes == 0] = 1
  return(sweep(U, 1, maxes, FUN="/"))
}
auc.concept <- function(concept)
{
  concept_gold_pairs = filter(gold_pairs, underlying_var_1 == concept) %>% 
    select(var_index_1, var_index_2)
  cos_sim_gold_pairs = apply(concept_gold_pairs, 1, retrieve.cos.sim) %>% na.omit()
  cos_sim_sample <- sample(cos_sim_matrix_melt$value, length(cos_sim_gold_pairs))
  yyi = c(cos_sim_gold_pairs, cos_sim_sample)
  Di = c(rep(1, length(cos_sim_gold_pairs)), rep(0, length(cos_sim_sample)))
  junk = ROC.Est.FUN(Di, yyi, yy0=0.5, fpr0=seq(0,1,0.01))
  return(junk[1])
}
auc.concept.cross <- function(concept)
{
  concept_gold_pairs = filter(gold_pairs, underlying_var_1 == concept)
  count0 = sum(concept_gold_pairs$split == 0)
  count1 = sum(concept_gold_pairs$split == 1)
  cos_sim_sample = c(sample(cos_sim_matrix_melt_0$value, count0), 
                     sample(cos_sim_matrix_melt_1$value, count1))
  yyi = c(concept_gold_pairs$cos_sim_gold_pairs, cos_sim_sample)
  Di = c(rep(1, nrow(concept_gold_pairs)), rep(0, length(cos_sim_sample)))
  junk = ROC.Est.FUN(Di, yyi, yy0=0.5, fpr0=seq(0,1,0.01))
  return(junk[1])
}

###########
#FILE PATHS
# CODER embeddings
file_chs = "coder_embeddings/desc-no-visit/chs_v7_coder_pp_embedding_desc-no-visit.csv"
file_mesa = "coder_embeddings/desc-no-visit/mesa_v13_coder_pp_embedding_desc-no-visit.csv"
file_whi = "coder_embeddings/desc-no-visit/whi_v12_coder_pp_embedding_desc-no-visit.csv"

# SapBERT embeddings
file_chs_sapbert = "coder_embeddings/sapbert_0609/sapbert_MEAN_chs_v7_v-desc-no-visit.csv"
file_mesa_sapbert = "coder_embeddings/sapbert_0609/sapbert_MEAN_mesa_v13_v-desc-no-visit.csv"
file_whi_sapbert = "coder_embeddings/sapbert_0609/sapbert_MEAN_whi_v12_v-desc-no-visit.csv"

# BioBERT embeddings
file_chs_biobert = "coder_embeddings/biobert_0609/biobert_MEAN_chs_v7_v-desc-no-visit.csv"
file_mesa_biobert = "coder_embeddings/biobert_0609/biobert_MEAN_mesa_v13_v-desc-no-visit.csv"
file_whi_biobert = "coder_embeddings/biobert_0609/biobert_MEAN_whi_v12_v-desc-no-visit.csv"


###########
# set splits
# hard concepts
concept_aucs = fread("auc_calc/intra_whi_concept_aucs_1007.csv")
gold_pairs_file <- "patient_data/gold_pairs/gold_pairs_intra_whi_0619.csv"
gold_pairs = fread(gold_pairs_file)
#colnames(concept_aucs) = c("underlying_var", "sapbert", "sonar")
underlying_vars = union(filter(concept_aucs, unsupervised < 0.9)$underlying_var, 
                        filter(concept_aucs, sap < 0.9)$underlying_var)
pairs_by_concept = table(gold_pairs$underlying_var_1) %>% data.frame()
nrow(pairs_by_concept)
table(pairs_by_concept$Freq == 1)

gold_pairs = gold_pairs %>% mutate(split = 0)

set.seed(1234)
onepair = filter(pairs_by_concept, Freq == 1) %>% filter(Var1 %in% underlying_vars)
onepair = onepair$Var1 %>% as.character() 
onepair = sample(onepair, size = round(length(onepair) / 2))

multpairs = filter(pairs_by_concept, Freq != 1) %>% filter(Var1 %in% underlying_vars)
multpairs = multpairs$Var1 %>% as.character()

change = c()
for(underlying_var in multpairs) {
  indices = which(gold_pairs$underlying_var_1 == underlying_var)
  change = c(change, sample(indices, size = round(length(indices) / 2)))
}
change = c(change, which(gold_pairs$underlying_var_1 %in% onepair))
gold_pairs$split[change] = 1

fwrite(gold_pairs, "patient_data/gold_pairs/gold_pairs_intra_whi_1017.csv")

# original
concept_aucs = fread("auc_calc/concept_aucs_mesa_1006.csv")
gold_pairs_file <- "patient_data/gold_pairs/gold_pairs_intra_mesa_0619.csv"
gold_pairs = fread(gold_pairs_file)
underlying_vars = unique(gold_pairs$underlying_var_1)
pairs_by_concept = table(gold_pairs$underlying_var_1) %>% data.frame()
nrow(pairs_by_concept)
table(pairs_by_concept$Freq == 1)

gold_pairs = gold_pairs %>% mutate(split = 0)

onepair = filter(pairs_by_concept, Freq == 1)
onepair = onepair$Var1 %>% as.character()
onepair = sample(onepair, size = round(length(onepair) / 2))

multpairs = filter(pairs_by_concept, Freq != 1)
multpairs = multpairs$Var1 %>% as.character()

change = c()
for(underlying_var in multpairs) {
  indices = which(gold_pairs$underlying_var_1 == underlying_var)
  change = c(change, sample(indices, size = round(length(indices) / 2)))
}
change = c(change, which(gold_pairs$underlying_var_1 %in% onepair))
gold_pairs$split[change] = 1

fwrite(gold_pairs, gold_pairs_file)

###########
# SET THESE
grp_means_clean_file = "data/var_means/inter_mesa-whi_quantiles_clean_grouped_0602.csv"
gold_pairs_file = "data/var_means/gold_pairs_inter_mesa-whi_0619.csv" #1010 for hard concepts, 0619 otherwise

###########
# SEMANTIC (same for intra and inter)
vars_wmeans = fread(grp_means_clean_file) %>% colnames()

# CHS
dat_chs_coder <- file_chs %>% fread()
var_key_chs <- data.frame(var_accession = dat_chs_coder[ , 1],
                          var_index = paste0("c", 1:nrow(dat_chs_coder)))
dat_chs_coder <- dat_chs_coder[ , -1] %>% data.matrix()
rownames(dat_chs_coder) <- var_key_chs$var_index
dat_chs_coder <- dat_chs_coder[rownames(dat_chs_coder) %in% vars_wmeans, ]
dat_chs_coder = rescale.max(dat_chs_coder)

dat_chs_sap <- file_chs_sapbert %>% fread()
var_key_chs <- data.frame(var_accession = dat_chs_sap[ , 1],
                          var_index = paste0("c", 1:nrow(dat_chs_sap)))
dat_chs_sap <- dat_chs_sap[ , -1] %>% data.matrix()
rownames(dat_chs_sap) <- var_key_chs$var_index
dat_chs_sap <- dat_chs_sap[rownames(dat_chs_sap) %in% vars_wmeans, ]
dat_chs_sap = rescale.max(dat_chs_sap)
# biobert
dat_chs_bio <- file_chs_biobert %>% fread()
var_key_chs <- data.frame(var_accession = dat_chs_bio[ , 1],
                          var_index = paste0("c", 1:nrow(dat_chs_bio)))
dat_chs_bio <- dat_chs_bio[ , -1] %>% data.matrix()
rownames(dat_chs_bio) <- var_key_chs$var_index
dat_chs_bio <- dat_chs_bio[rownames(dat_chs_bio) %in% vars_wmeans, ]

# MESA
dat_mesa_coder <- file_mesa %>% fread()
var_key_mesa <- data.frame(var_accession = dat_mesa_coder[ , 1],
                           var_index = paste0("m", 1:nrow(dat_mesa_coder)))
dat_mesa_coder <- dat_mesa_coder[ , -1] %>% data.matrix()
rownames(dat_mesa_coder) <- var_key_mesa$var_index
dat_mesa_coder <- dat_mesa_coder[rownames(dat_mesa_coder) %in% vars_wmeans, ]
dat_mesa_coder = rescale.max(dat_mesa_coder)

dat_mesa_sap <- file_mesa_sapbert %>% fread()
var_key_mesa <- data.frame(var_accession = dat_mesa_sap[ , 1],
                           var_index = paste0("m", 1:nrow(dat_mesa_sap)))
dat_mesa_sap <- dat_mesa_sap[ , -1] %>% data.matrix()
rownames(dat_mesa_sap) <- var_key_mesa$var_index
dat_mesa_sap <- dat_mesa_sap[rownames(dat_mesa_sap) %in% vars_wmeans, ]
dat_mesa_sap = rescale.max(dat_mesa_sap)
# biobert
dat_mesa_bio <- file_mesa_biobert %>% fread()
var_key_mesa <- data.frame(var_accession = dat_mesa_bio[ , 1],
                           var_index = paste0("m", 1:nrow(dat_mesa_bio)))
dat_mesa_bio <- dat_mesa_bio[ , -1] %>% data.matrix()
rownames(dat_mesa_bio) <- var_key_mesa$var_index
dat_mesa_bio <- dat_mesa_bio[rownames(dat_mesa_bio) %in% vars_wmeans, ]

# WHI
dat_whi_coder <- file_whi %>% fread()
var_key_whi <- data.frame(var_accession = dat_whi_coder[ , 1],
                          var_index = paste0("w", 1:nrow(dat_whi_coder)))
dat_whi_coder <- dat_whi_coder[ , -1] %>% data.matrix()
rownames(dat_whi_coder) <- var_key_whi$var_index
dat_whi_coder <- dat_whi_coder[rownames(dat_whi_coder) %in% vars_wmeans, ]
dat_whi_coder = rescale.max(dat_whi_coder)

dat_whi_sap <- file_whi_sapbert %>% fread()
var_key_whi <- data.frame(var_accession = dat_whi_sap[ , 1],
                          var_index = paste0("w", 1:nrow(dat_whi_sap)))
dat_whi_sap <- dat_whi_sap[ , -1] %>% data.matrix()
rownames(dat_whi_sap) <- var_key_whi$var_index
dat_whi_sap <- dat_whi_sap[rownames(dat_whi_sap) %in% vars_wmeans, ]
dat_whi_sap = rescale.max(dat_whi_sap)
# biobert
dat_whi_bio <- file_whi_biobert %>% fread()
var_key_whi <- data.frame(var_accession = dat_whi_bio[ , 1],
                          var_index = paste0("w", 1:nrow(dat_whi_bio)))
dat_whi_bio <- dat_whi_bio[ , -1] %>% data.matrix()
rownames(dat_whi_bio) <- var_key_whi$var_index
dat_whi_bio <- dat_whi_bio[rownames(dat_whi_bio) %in% vars_wmeans, ]


###########
# INTRA
# DISTRIBUTION
chs_means <- fread(grp_means_clean_file) %>% t()
chs_means = rescale.max(chs_means)
chs_means <- chs_means[match(rownames(dat_chs_coder), rownames(chs_means)), ]
chs = cbind(dat_chs_coder, dat_chs_sap, chs_means)

mesa_means <- fread(grp_means_clean_file) %>% t()
mesa_means = rescale.max(mesa_means)
mesa_means <- mesa_means[match(rownames(dat_mesa_coder), rownames(mesa_means)), ]
mesa = cbind(dat_mesa_coder, dat_mesa_sap, mesa_means)

whi_means <- fread(grp_means_clean_file) %>% t()
whi_means = rescale.max(whi_means)
whi_means <- whi_means[match(rownames(dat_whi_coder), rownames(whi_means)), ]
whi = cbind(dat_whi_coder, dat_whi_sap, whi_means)

# COS SIM MATRICES for xgboost
cos_sim_matrix_sonar = cos.sim.matrix(chs)
cos_sim_matrix_sap = cos.sim.matrix(dat_chs_sap)
cos_sim_matrix_coder = cos.sim.matrix(dat_chs_coder)
cos_sim_matrix_dist = cos.sim.matrix(chs_means)

cos_sim_matrix_sonar = cos.sim.matrix(mesa)
cos_sim_matrix_sap = cos.sim.matrix(dat_mesa_sap)
cos_sim_matrix_coder = cos.sim.matrix(dat_mesa_coder)
cos_sim_matrix_dist = cos.sim.matrix(mesa_means)

cos_sim_matrix_sonar = cos.sim.matrix(whi)
cos_sim_matrix_sap = cos.sim.matrix(dat_whi_sap)
cos_sim_matrix_coder = cos.sim.matrix(dat_whi_coder)
cos_sim_matrix_dist = cos.sim.matrix(whi_means)

# LOSS FUNCTION
gold_pairs = fread(gold_pairs_file)
underlying_vars = unique(gold_pairs$underlying_var_1)
gold_pairs_0 = filter(gold_pairs, split == 0)
gold_pairs_1 = filter(gold_pairs, split == 1)

intra = whi
# trained on split = 0
gold_pairs_2 = data.table(other=gold_pairs_0$var_index_1,
                          loinc=gold_pairs_0$var_index_2,
                          ans=rep(1, nrow(gold_pairs_0)))
set.seed(123)
sample = data.table(other=sample(rownames(intra), nrow(gold_pairs_2)), 
                    loinc=sample(rownames(intra), nrow(gold_pairs_2)),
                    ans=rep(0, nrow(gold_pairs_2)))
dict_label = rbind(gold_pairs_2, sample)
start.time <- proc.time() # start time
intra_new = get_supervised(Other=intra, LOINC=intra, dict_label=dict_label, coef = c(2,50,0.1), 
                           stepsize = 0.001, maxstep = 200, epsilon = 0.01)[[1]]
proc.time() - start.time
proc.time()# end time

intra_new = fread("supervised_embeddings/intra_whi_concat_split0_0910.csv")
intra_new = as.matrix(intra_new)
rownames(intra_new) = rownames(dat_whi_coder)

cos_sim_matrix = cos.sim.matrix(intra_new)
cos_sim_matrix_0 = cos_sim_matrix
cos_sim_matrix_melt_0 <- melt.data.table(data.table(cos_sim_matrix)) %>% na.omit()
cos_sim_gold_pairs_0 = apply(select(gold_pairs_1, var_index_1, var_index_2), 1, retrieve.cos.sim)

# trained on split = 1
gold_pairs_2 = data.table(other=gold_pairs_1$var_index_1,
                          loinc=gold_pairs_1$var_index_2,
                          ans=rep(1, nrow(gold_pairs_1)))
set.seed(12)
sample = data.table(other=sample(rownames(intra), nrow(gold_pairs_2)), 
                    loinc=sample(rownames(intra), nrow(gold_pairs_2)),
                    ans=rep(0, nrow(gold_pairs_2)))
dict_label = rbind(gold_pairs_2, sample)
start.time <- proc.time() # start time
intra_new = get_supervised(Other=intra, LOINC=intra, dict_label=dict_label, coef = c(2,50,0.1), 
                           stepsize = 0.001, maxstep = 200, epsilon = 0.01)[[1]] # try this
proc.time() - start.time
proc.time()# end time

intra_new = fread("supervised_embeddings/intra_whi_concat_split1_0910.csv")
intra_new = as.matrix(intra_new)
rownames(intra_new) = rownames(dat_whi_coder)

cos_sim_matrix = cos.sim.matrix(intra_new)
cos_sim_matrix_1 = cos_sim_matrix
cos_sim_matrix_melt_1 <- melt.data.table(data.table(cos_sim_matrix)) %>% na.omit()
cos_sim_gold_pairs_1 = apply(select(gold_pairs_0, var_index_1, var_index_2), 1, retrieve.cos.sim)


###########
# INTER
# DISTRIBUTION
means <- fread(grp_means_clean_file) %>% t()
means = rescale.max(means)

chs_vars <- grepl("c", rownames(means))
chs_means <- means[chs_vars, ]
chs_means <- chs_means[match(rownames(dat_chs_coder), rownames(chs_means)), ]
chs = cbind(dat_chs_coder, dat_chs_sap, chs_means)

mesa_vars <- grepl("m", rownames(means))
mesa_means <- means[mesa_vars, ]
mesa_means <- mesa_means[match(rownames(dat_mesa_coder), rownames(mesa_means)), ]
mesa = cbind(dat_mesa_coder, dat_mesa_sap, mesa_means)

whi_vars <- grepl("w", rownames(means))
whi_means <- means[whi_vars, ]
whi_means <- whi_means[match(rownames(dat_whi_coder), rownames(whi_means)), ]
whi = cbind(dat_whi_coder, dat_whi_sap, whi_means)

cos_sim_matrix_sonar = cos.sim.two.matrices(chs, mesa)
cos_sim_matrix_sap = cos.sim.two.matrices(dat_chs_sap, dat_mesa_sap)
cos_sim_matrix_coder = cos.sim.two.matrices(dat_chs_coder, dat_mesa_coder)
cos_sim_matrix_dist = cos.sim.two.matrices(chs_means, mesa_means)

cos_sim_matrix_sonar = cos.sim.two.matrices(chs, whi)
cos_sim_matrix_sap = cos.sim.two.matrices(dat_chs_sap, dat_whi_sap)
cos_sim_matrix_coder = cos.sim.two.matrices(dat_chs_coder, dat_whi_coder)
cos_sim_matrix_dist = cos.sim.two.matrices(chs_means, whi_means)

cos_sim_matrix_sonar = cos.sim.two.matrices(mesa, whi)
cos_sim_matrix_sap = cos.sim.two.matrices(dat_mesa_sap, dat_whi_sap)
cos_sim_matrix_coder = cos.sim.two.matrices(dat_mesa_coder, dat_whi_coder)
cos_sim_matrix_dist = cos.sim.two.matrices(mesa_means, whi_means)

# LOSS FUNCTION
gold_pairs = fread(gold_pairs_file)
underlying_vars = unique(gold_pairs$underlying_var_1)
gold_pairs_0 = filter(gold_pairs, split == 0)
gold_pairs_1 = filter(gold_pairs, split == 1)

inter1 = mesa
inter2 = whi
# trained on split = 0
gold_pairs_2 = data.table(other=gold_pairs_0$var_index_1,
                          loinc=gold_pairs_0$var_index_2,
                          ans=rep(1, nrow(gold_pairs_0)))
set.seed(123)
sample = data.table(other=sample(rownames(inter1), nrow(gold_pairs_2)), 
                    loinc=sample(rownames(inter2), nrow(gold_pairs_2)),
                    ans=rep(0, nrow(gold_pairs_2)))
dict_label = rbind(gold_pairs_2, sample)

start.time <- proc.time() # start time
inter1_new = get_supervised(Other=inter1, LOINC=inter2, dict_label=dict_label, coef = c(2,50,0.1), 
                            stepsize = 0.001, maxstep = 200)[[1]]
proc.time() - start.time
proc.time()# end time

cos_sim_matrix = cos.sim.two.matrices(inter1_new, inter2)
cos_sim_matrix_0 = cos_sim_matrix
cos_sim_matrix_melt_0 <- melt.data.table(data.table(cos_sim_matrix)) %>% na.omit()
cos_sim_gold_pairs_0 = apply(select(gold_pairs_1, var_index_1, var_index_2), 1, retrieve.cos.sim)

# trained on split = 1
gold_pairs_2 = data.table(other=gold_pairs_1$var_index_1,
                          loinc=gold_pairs_1$var_index_2,
                          ans=rep(1, nrow(gold_pairs_1)))
set.seed(12)
sample = data.table(other=sample(rownames(inter1), nrow(gold_pairs_2)), 
                    loinc=sample(rownames(inter2), nrow(gold_pairs_2)),
                    ans=rep(0, nrow(gold_pairs_2)))
dict_label = rbind(gold_pairs_2, sample)

start.time <- proc.time() # start time
inter1_new = get_supervised(Other=inter1, LOINC=inter2, dict_label=dict_label, coef = c(2,50,0.1), 
                            stepsize = 0.001, maxstep = 200)[[1]] # try this
proc.time() - start.time
proc.time()# end time

cos_sim_matrix = cos.sim.two.matrices(inter1_new, inter2)
cos_sim_matrix_1 = cos_sim_matrix
cos_sim_matrix_melt_1 <- melt.data.table(data.table(cos_sim_matrix)) %>% na.omit()
cos_sim_gold_pairs_1 = apply(select(gold_pairs_0, var_index_1, var_index_2), 1, retrieve.cos.sim)


###########
# AUC (WITHOUT CROSS-VALIDATION)
gold_pairs = fread(gold_pairs_file)
underlying_vars = unique(gold_pairs$underlying_var_1)
cos_sim_matrix_melt <- melt.data.table(data.table(cos_sim_matrix)) %>% na.omit()
aucs = c()
for(seed in c(12, 123, 1234)) {
  set.seed(seed)
  for (concept in underlying_vars) {
    aucs = c(aucs, auc.concept(concept))
  }
}
mean(aucs)

aucs = matrix(aucs, ncol=3)
aucs = apply(aucs, 1, mean)
concept_aucs = cbind(concept_aucs, aucs)
colnames(concept_aucs) = c("underlying_var", "sap")
colnames(concept_aucs) = c("underlying_var", "supervised", "unsupervised", "dist", "coder_sap", "coder", "sap", "bio")

# AUC (WITH CROSS-VALIDATION)
cos_sim_gold_pairs = c(cos_sim_gold_pairs_0, cos_sim_gold_pairs_1)
gold_pairs = rbind(gold_pairs_1, gold_pairs_0) %>% cbind(cos_sim_gold_pairs)
cos_sim_matrix = Reduce("+", list(cos_sim_matrix_0, cos_sim_matrix_1)) / 2
for(i in 1:nrow(gold_pairs)) {
  cos_sim_matrix[gold_pairs[i,]$var_index_1, gold_pairs[i,]$var_index_2] = 
    gold_pairs[i,]$cos_sim_gold_pairs
  # for intra
  cos_sim_matrix[gold_pairs[i,]$var_index_2, gold_pairs[i,]$var_index_1] = 
    gold_pairs[i,]$cos_sim_gold_pairs
}

aucs = c()
for(seed in c(1000, 2000, 3000)) {
  set.seed(seed)
  for (concept in underlying_vars) {
    aucs = c(aucs, auc.concept.cross(concept))
  }
}
mean(aucs)

###########
# GET SENSITIVITY

# WITHOUT CROSS VALIDATION
gold_pairs = fread(gold_pairs_file)
gold_vars = unique(c(gold_pairs$var_index_1, gold_pairs$var_index_2))

# WITH CROSS VALIDATION
cos_sim_matrix = Reduce("+", list(cos_sim_matrix_0, cos_sim_matrix_1)) / 2
for(i in 1:nrow(gold_pairs)) {
  cos_sim_matrix[gold_pairs[i,]$var_index_1, gold_pairs[i,]$var_index_2] = 
    gold_pairs[i,]$cos_sim_gold_pairs
  # for intra
  cos_sim_matrix[gold_pairs[i,]$var_index_2, gold_pairs[i,]$var_index_1] = 
    gold_pairs[i,]$cos_sim_gold_pairs
}
gold_vars = unique(c(gold_pairs$var_index_1, gold_pairs$var_index_2))

# INTRA
diag(cos_sim_matrix) = 0
cos_sim_gold_pairs = cos_sim_matrix[, gold_vars]

top_50 = apply(cos_sim_gold_pairs, 2, kit::topn, 50)
key = data.frame(var_index = rownames(cos_sim_matrix), topn_index = 1:nrow(cos_sim_matrix))

top_50 %>% write.csv("data/var_means/top50_intra_whi_biobert_0628.csv", row.names = FALSE)
key %>% write.csv("data/var_means/top50_KEY_intra_whi_0712.csv", row.names = FALSE)

# INTER
cos_sim_gold_pairs = t(cos_sim_matrix[unique(gold_pairs$var_index_1), ])
top_50 = apply(cos_sim_gold_pairs, 2, kit::topn, 50)

cos_sim_gold_pairs = cos_sim_matrix[, unique(gold_pairs$var_index_2)]
top_50 = cbind(top_50, apply(cos_sim_gold_pairs, 2, kit::topn, 50))

key = data.frame(var_index = c(rownames(cos_sim_matrix), colnames(cos_sim_matrix)), 
                 topn_index = c(1:nrow(cos_sim_matrix), 1:ncol(cos_sim_matrix)))

top_50 %>% write.csv("data/var_means/top50_inter_mesa-whi_biobert_0628.csv", row.names = FALSE)
key %>% write.csv("data/var_means/top50_KEY_inter_mesa-whi4_0717.csv", row.names = FALSE)

###########
# SENSITIVITY CALC (INTRA)
top_50 = read.csv("patient_data/topk/top50_intra_mesa_xgb_1029.csv")
top_50_key = read.csv("patient_data/topk/top50_KEY_intra_mesa_1029.csv")
# "patient_data/topk/top50_intra_whi_lossfxn_0801.csv"
# "patient_data/topk/top50_intra_whi_lossfxn_0720.csv"
# "patient_data/topk/top50_intra_whi_dist_0628.csv"
# "patient_data/topk/top50_intra_whi_biobert_0612.csv"
# "patient_data/topk/top50_KEY_intra_whi_0628.csv"

gold_pairs = fread("patient_data/gold_pairs/gold_pairs_intra_mesa_0619.csv")
gold_pairs = gold_pairs %>% select(underlying_var_1, var_index_1, var_index_2)
gold_pairs_2 = select(gold_pairs, underlying_var_1, var_index_2, var_index_1)
colnames(gold_pairs_2) = c("underlying_var_1", "var_index_1", "var_index_2")
gold_pairs = rbind(gold_pairs, gold_pairs_2)
gold_pairs = gold_pairs %>% left_join(top_50_key, by=c("var_index_2" = "var_index"))

gold_vars = cbind(c(gold_pairs$underlying_var_1, gold_pairs$underlying_var_1), 
                  c(gold_pairs$var_index_1, gold_pairs$var_index_2)) %>% data.table()
colnames(gold_vars) = c("underlying_var", "var_index")

gold_vars = gold_vars %>% distinct(var_index, .keep_all = TRUE)
gold_vars = data.table(var_index = colnames(top_50)) %>% left_join(gold_vars, by = "var_index")

# this one
topk = c(1, 3, 5, 10, 20)
junk = data.frame(underlying_var = c(), percent_not0 = c(), k = c())
for(k in topk) {
  percent_matches = c()
  for(var in colnames(top_50)) {
    var_top = top_50[, var]
    var_matches = filter(gold_pairs, var_index_1 == var)
    var_top = var_top[1:k]
    percent_var_matches = sum(var_matches$topn_index %in% var_top) / nrow(var_matches)
    percent_matches = c(percent_matches, percent_var_matches)
  }
  percent_matches = cbind(gold_vars, percent_matches) %>%
    mutate(not0 = (percent_matches != 0)) %>%
    group_by(underlying_var) %>%
    summarise(percent_not0 = mean(not0)) %>%
    mutate(top = k)
  colnames(percent_matches) = c("underlying_var", "percent_not0", "k")
  mean(percent_matches$percent_not0)
  junk = rbind(junk, percent_matches)
}

junk %>% group_by(k) %>% summarise(avg_percent_not0 = mean(percent_not0))
junk = mutate(junk, comparison = "whi", comparison_cat = "intra")

#all = junk
all = rbind(all, junk)


# SENSITIVITY CALC (INTER)
top_50 = read.csv("patient_data/topk/top50_inter_chs-mesa_xgb_1029.csv")
top_50_key = read.csv("patient_data/topk/top50_KEY_inter_chs-mesa_0717.csv")
# "patient_data/topk/top50_inter_mesa-whi_lossfxn3_0717.csv"
# "patient_data/topk/top50_inter_mesa-whi_dist_0628.csv"

gold_pairs = fread("patient_data/gold_pairs/gold_pairs_inter_chs-mesa_0619.csv")
gold_pairs = gold_pairs %>% select(underlying_var_1, var_index_1, var_index_2)

# to reverse direction
colnames(gold_pairs) = c("underlying_var_1", "var_index_2", "var_index_1")

gold_pairs = gold_pairs %>% left_join(top_50_key, by=c("var_index_2" = "var_index"))

gold_vars = select(gold_pairs, underlying_var_1, var_index_1) %>% 
  distinct(var_index_1, .keep_all = TRUE) %>%
  filter(var_index_1 %in% colnames(top_50))

gold_pairs = gold_pairs %>% filter(var_index_1 %in% colnames(top_50)) %>%
  filter(var_index_2 %in% colnames(top_50))

# this one
topk = c(1, 3, 5, 10, 20)
junk = data.frame(underlying_var = c(), percent_not0 = c(), k = c())
for(k in topk) {
  percent_matches = c()
  for(var in unique(gold_pairs$var_index_1)) {
    var_top = top_50[, var]
    var_matches = filter(gold_pairs, var_index_1 == var)
    var_top = var_top[1:k]
    percent_var_matches = sum(var_matches$topn_index %in% var_top) / nrow(var_matches)
    percent_matches = c(percent_matches, percent_var_matches)
  }
  percent_matches = cbind(gold_vars, percent_matches) %>%
    mutate(not0 = (percent_matches != 0)) %>%
    group_by(underlying_var_1) %>%
    summarise(percent_not0 = mean(not0)) %>%
    mutate(top = k)
  colnames(percent_matches) = c("underlying_var", "percent_not0", "k")
  mean(percent_matches$percent_not0)
  junk = rbind(junk, percent_matches)
}
junk %>% group_by(k) %>% summarise(avg_percent_not0 = mean(percent_not0))
junk = mutate(junk, comparison = "mesa-chs", comparison_cat = "inter")

#all = junk
all = rbind(all, junk)

all %>% write.csv("patient_data/topk/all_xgb_1029.csv", row.names = F)

old = options(pillar.sigfig = 4)

### THIS
junk = all %>% group_by(k, underlying_var) %>% summarise(avg_percent_not0 = mean(percent_not0))
junk = junk %>% group_by(k) %>% summarise(avg_percent_not0 = mean(avg_percent_not0))
junk = cbind("all", junk)
colnames(junk) = c("comparison", "k", "avg_percent_not0")
save = junk

junk = all %>% group_by(comparison_cat, k, underlying_var) %>% 
  summarise(avg_percent_not0 = mean(percent_not0))
junk = junk %>% group_by(comparison_cat, k) %>% summarise(avg_percent_not0 = mean(avg_percent_not0))
colnames(junk) = c("comparison", "k", "avg_percent_not0")
save = rbind(save, junk)

junk = all %>% group_by(comparison, k, underlying_var) %>% 
  summarise(avg_percent_not0 = mean(percent_not0))
junk = junk %>% group_by(comparison, k) %>% summarise(avg_percent_not0 = mean(avg_percent_not0))
save = rbind(save, junk)
colnames(save) = c("comparison", "k", "xgb")

#all_sum = save
all_sum = left_join(all_sum, save, by = c("comparison", "k"))

all_sum %>% write.csv("patient_data/topk/all_summary_1029.csv", row.names = F)

t(select(all_sum, comparison, k, dist, biobert, coder, sap, "coder-sap", "coder-sap-max", 
         untrained, combined, xgb)) %>% 
  write.csv("patient_data/topk/all_summary_transpose_1029.csv", row.names = T)



all_copy = all

all %>% group_by(comparison_cat, k) %>% summarise(avg_percent_not0 = mean(percent_not0))

all %>% group_by(k, underlying_var) %>% summarise(avg_percent_not0 = mean(percent_not0))


junk = rbind(junk, percent_matches)
junk = percent_matches

junk = junk %>% group_by(underlying_var) %>% summarise(avg_percent_not0 = mean(percent_not0))

save = left_join(save, junk, by = "underlying_var")
colnames(save) = c("underlying_var", "top1", "top3", "top5", "top10", "top20",
                   "sap_top1", "sap_top3", "sap_top5", "sap_top10", "sap_top20",
                   "coder_top1", "coder_top3", "coder_top5", "coder_top10", "coder_top20",
                   "coder_sap_top1","coder_sap_top3","coder_sap_top5","coder_sap_top10","coder_sap_top20",
                   "dist_top1", "dist_top3", "dist_top5", "dist_top10", "dist_top20")

######
# Pool all gold-standard pairs
gold_pairs_intra = rbind(read.csv("patient_data/gold_pairs/gold_pairs_intra_chs_0619.csv"),
                         read.csv("patient_data/gold_pairs/gold_pairs_intra_mesa_0619.csv"),
                         read.csv("patient_data/gold_pairs/gold_pairs_intra_whi_0619.csv")) %>% 
  select(underlying_var_1, var_index_1, var_index_2)
gold_pairs_inter = rbind(read.csv("patient_data/gold_pairs/gold_pairs_inter_chs-mesa_0619.csv"),
                         read.csv("patient_data/gold_pairs/gold_pairs_inter_chs-whi_0619.csv"),
                         read.csv("patient_data/gold_pairs/gold_pairs_inter_mesa-whi_0619.csv")) %>% 
  select(underlying_var_1, var_index_1, var_index_2)

gold_pairs = rbind(gold_pairs_intra, gold_pairs_inter)
