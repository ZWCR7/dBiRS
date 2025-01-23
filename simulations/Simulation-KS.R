 

library(BiRS)
library(QSCAN)
library(genio)
library(Matrix)
library(MASS)
library(foreach)
library(seqminer)
library(doParallel)
library(KnockoffScreen)

BIMList = read_bim('Gene_Sample/genotype.bim')

chr = 1
BP = BIMList$pos; index_trunc = max(which(BP < 5000*1000))

BP = BP[1:index_trunc]
RSID = BIMList$id[1:index_trunc]

posmin = min(BP); posmax = max(BP)

pos.tag = seq(posmin, posmax, by = 5000/2)
window.bed = cbind(chr, pos.tag, pos.tag + 5000 - 3)
window.bed = window.bed[order(as.numeric(window.bed[, 2])),]

ind.start = 1; Block.bed = NULL; Block.ind = NULL
pos.end = 0

while(pos.end < posmax)
{
  pos.start = BP[ind.start]
  pos.trunc = pos.start + 50*1000
  
  ind.end = max(which(BP <= pos.trunc))
  pos.end = BP[ind.end]
  
  Block.bed = rbind(Block.bed, c(chr, pos.start, pos.end))
  Block.ind = rbind(Block.ind, c(chr, ind.start, ind.end))
  
  print(ind.end - ind.start)
  ind.start = ind.end + 1
}

set.seed(1024)

num_region = 4
locuind = c(15, 36, 60, 84)

p = length(RSID)
FAMList = read_fam('Gene_Sample/genotype.fam')

n = 10000
select_sample = sort(sample(1:length(FAMList$fam), size = n, replace = F)); FinalID = rep(0, n)
for (i in 1:n)
{
  FinalID[i] = paste0(FAMList$fam[select_sample[i]], '_', FAMList$id[select_sample[i]])
}

X1 = rnorm(n, 0, 1); X2 = sample(c(0, 1), size = n, replace = T)
X = cbind(X1, X2)

gamma = c(0.5, 0.5)

GenoFile = 'Gene_Sample/genotype.vcf.gz'
GenoPrefix = 'Gene_Sample/genotype'

source('Impute_Func.R')

beta = rep(0, p)
set.seed(2)
Yglobal = rep(0, n)
index_causal = c()
for (i in 1:num_region)
{
  Gi = -1*(readPlinkToMatrixByIndex(plinkFilePrefix = GenoPrefix, sampleIndex = select_sample, markerIndex = Block.ind[locuind[i], 2]:Block.ind[locuind[i], 3]) - 2)
  Gi[Gi < 0 | Gi > 2] = NA
  
  if (sum(is.na(Gi)) > 0)
  {
    Gi = apply(Gi, 2, Impute_v, impute.method = 'bestguess')
  }
  
  MAFi = colMeans(Gi)/2; #EleVari = apply(Gi, 2, var)
  #print(sum(2*MAFi < 0.001))
  MAFi[which(MAFi == 0)] = 0.0001
  colBP = as.numeric(gsub("^.*\\:","",colnames(Gi)))
  
  betai = rep(0, length(colBP))
  
  causal.posstart = sample(colBP, 1); causal.posend = causal.posstart + 5*1000
  causal.indstart = which(colBP == causal.posstart); causal.indend = max(which(colBP < causal.posend))
  
  signalpos = causal.indstart:causal.indend
  
  signalpos1 = sample(signalpos, size = ceiling(0.1*length(signalpos)), replace = F)
  signalpos2 = setdiff(signalpos, signalpos1)
  signsel = sample(c(-1, 1), size = length(signalpos1), prob = c(0.5, 0.5), replace = T)
  
  print(MAFi[signalpos1])
  
  betai[signalpos1] = abs(log10(MAFi[signalpos1]))*signsel
  
  betai[signalpos2] = 0
  
  origin.startind = which(BP == min(colBP)); origin.endind = which(BP == max(colBP))
  index_causal = c(index_causal, origin.startind + signalpos - 1)
  
  print(length(signalpos))
  
  Yglobal = Yglobal + Gi %*% betai
  beta[origin.startind:origin.endind] = betai
}

rm(Gi); gc()

cv = c(0.12, 0.15)
for (setc in 1:length(cv))
{
  c = cv[setc]
  Ycore = X %*% gamma + c*Yglobal
  
  nsimu = 100
  
  Simu_dBiRS = function(s)
  {
    set.seed(s)
    
    epsilon = rnorm(n, 0, 1)
    
    Y = Ycore + epsilon
    phenotype = as.vector(Y)
    
    for (part in 1:25)
    {
      start_ind = 200*1000*(part - 1) + 1; end_ind = 200*1000*part
      posmin = min(BP[which(BP >= start_ind)]); posmax = max(BP[which(BP <= end_ind)])
      
      BP_part = posmin:posmax
      
      pos.tag = seq(posmin, posmax, by = 5000/2)
      window.bed = cbind(chr, pos.tag, pos.tag + 5000 - 3)
      window.bed = window.bed[order(as.numeric(window.bed[, 2])),]
      
      result.prelim = KS.prelim(phenotype, X = X, id = FinalID, out_type = "C")
      
      result.fit = KS.chr(result.prelim = result.prelim, seq.filename = GenoFile, window.bed = window.bed, region.pos = NULL,
                          tested.pos = BP_part, excluded.pos = NULL, M = 5, thres.single = 0.01, thres.ultrarare = 1, thres.missing = 0.05,
                          midout.dir = NULL, temp.dir = NULL, jobtitle = NULL, Gsub.id = NULL, impute.method = "bestguess", bigmemory = T,
                          leveraging = T, LD.filter = 0.75)
      
      save(list = c('result.fit'), file = paste0('result_ks_c1/KS_c = ', c, 's = ', s, 'part', part, '.RData'))
      
    }
  }
  
  print('KS Detection-----------------Start')
  aa = Sys.time()
  cl = makeCluster(100)
  registerDoParallel(cl)
  
  RES_BBiRS = foreach(s = 1:nsimu, .packages = c("BiRS", "MASS", "QSCAN", "Matrix", "seqminer", "genio", "KnockoffScreen")) %dopar% Simu_dBiRS(s)
  
  stopImplicitCluster()
  stopCluster(cl)
  bb = Sys.time()
  print(paste0('KS Detection', 'c' = c, '---------------End, Time = ', bb - aa))
}
########################################################################################################################################################################################
rm(list = ls())
library(BiRS)
library(QSCAN)
library(genio)
library(Matrix)
library(MASS)
library(foreach)
library(seqminer)
library(doParallel)
library(KnockoffScreen)

BIMList = read_bim('Gene_Sample/genotype.bim')

chr = 1
BP = BIMList$pos; index_trunc = max(which(BP < 5000*1000))

BP = BP[1:index_trunc]
RSID = BIMList$id[1:index_trunc]

posmin = min(BP); posmax = max(BP)

pos.tag = seq(posmin, posmax, by = 5000/2)
window.bed = cbind(chr, pos.tag, pos.tag + 5000 - 3)
window.bed = window.bed[order(as.numeric(window.bed[, 2])),]

ind.start = 1; Block.bed = NULL; Block.ind = NULL
pos.end = 0

while(pos.end < posmax)
{
  pos.start = BP[ind.start]
  pos.trunc = pos.start + 50*1000
  
  ind.end = max(which(BP <= pos.trunc))
  pos.end = BP[ind.end]
  
  Block.bed = rbind(Block.bed, c(chr, pos.start, pos.end))
  Block.ind = rbind(Block.ind, c(chr, ind.start, ind.end))
  
  print(ind.end - ind.start)
  ind.start = ind.end + 1
}

set.seed(1024)

num_region = 4
locuind = c(15, 36, 60, 84)

p = length(RSID)
FAMList = read_fam('Gene_Sample/genotype.fam')

n = 10000
select_sample = sort(sample(1:length(FAMList$fam), size = n, replace = F)); FinalID = rep(0, n)
for (i in 1:n)
{
  FinalID[i] = paste0(FAMList$fam[select_sample[i]], '_', FAMList$id[select_sample[i]])
}

X1 = rnorm(n, 0, 1); X2 = sample(c(0, 1), size = n, replace = T)
X = cbind(X1, X2)

gamma = c(0.5, 0.5)

GenoFile = 'Gene_Sample/genotype.vcf.gz'
GenoPrefix = 'Gene_Sample/genotype'

source('Impute_Func.R')

beta = rep(0, p)
set.seed(2)
Yglobal = rep(0, n)
index_causal = c()
for (i in 1:num_region)
{
  Gi = -1*(readPlinkToMatrixByIndex(plinkFilePrefix = GenoPrefix, sampleIndex = select_sample, markerIndex = Block.ind[locuind[i], 2]:Block.ind[locuind[i], 3]) - 2)
  Gi[Gi < 0 | Gi > 2] = NA
  
  if (sum(is.na(Gi)) > 0)
  {
    Gi = apply(Gi, 2, Impute_v, impute.method = 'bestguess')
  }
  
  MAFi = colMeans(Gi)/2; #EleVari = apply(Gi, 2, var)
  #print(sum(2*MAFi < 0.001))
  MAFi[which(MAFi == 0)] = 0.0001
  colBP = as.numeric(gsub("^.*\\:","",colnames(Gi)))
  
  betai = rep(0, length(colBP))
  
  causal.posstart = sample(colBP, 1); causal.posend = causal.posstart + 5*1000
  causal.indstart = which(colBP == causal.posstart); causal.indend = max(which(colBP < causal.posend))
  
  signalpos = causal.indstart:causal.indend
  
  signalpos1 = sample(signalpos, size = ceiling(0.1*length(signalpos)), replace = F)
  signalpos2 = setdiff(signalpos, signalpos1)
  signsel = sample(c(-1, 1), size = length(signalpos1), prob = c(0.5, 0.5), replace = T)
  
  print(MAFi[signalpos1])
  
  betai[signalpos1] = abs(log10(MAFi[signalpos1]))*signsel
  
  betai[signalpos2] = 0
  
  origin.startind = which(BP == min(colBP)); origin.endind = which(BP == max(colBP))
  index_causal = c(index_causal, origin.startind + signalpos - 1)
  
  print(length(signalpos))
  
  Yglobal = Yglobal + Gi %*% betai
  beta[origin.startind:origin.endind] = betai
}

rm(Gi); gc()

inv_logi = function(x) exp(x)/(1 + exp(x))

cv = c(0.25, 0.30)
for (setc in 1:length(cv))
{
  c = cv[setc]
  Ycore = X %*% gamma + c*Yglobal
  
  nsimu = 100
  
  Simu_dBiRS = function(s)
  {
    set.seed(s)
    
    eta = Ycore
    pi = inv_logi(eta)
    Y = rbinom(n, size = 1, prob = pi)
    phenotype = as.vector(Y)
    
    for (part in 1:25)
    {
      start_ind = 200*1000*(part - 1) + 1; end_ind = 200*1000*part
      posmin = min(BP[which(BP >= start_ind)]); posmax = max(BP[which(BP <= end_ind)])
      
      BP_part = posmin:posmax
      
      pos.tag = seq(posmin, posmax, by = 5000/2)
      window.bed = cbind(chr, pos.tag, pos.tag + 5000 - 3)
      window.bed = window.bed[order(as.numeric(window.bed[, 2])),]
      
      result.prelim = KS.prelim(phenotype, X = X, id = FinalID, out_type = "D")
      
      result.fit = KS.chr(result.prelim = result.prelim, seq.filename = GenoFile, window.bed = window.bed, region.pos = NULL,
                          tested.pos = BP_part, excluded.pos = NULL, M = 5, thres.single = 0.01, thres.ultrarare = 1, thres.missing = 0.05,
                          midout.dir = NULL, temp.dir = NULL, jobtitle = NULL, Gsub.id = NULL, impute.method = "bestguess", bigmemory = T,
                          leveraging = T, LD.filter = 0.75)
      
      save(list = c('result.fit'), file = paste0('result_ks_b1/KS_c = ', c, 's = ', s, 'part', part, '.RData'))
      
    }
  }
  
  print('KS Detection-----------------Start')
  aa = Sys.time()
  cl = makeCluster(100)
  registerDoParallel(cl)
  
  RES_BBiRS = foreach(s = 1:nsimu, .packages = c("BiRS", "MASS", "QSCAN", "Matrix", "seqminer", "genio", "KnockoffScreen")) %dopar% Simu_dBiRS(s)
  
  stopImplicitCluster()
  stopCluster(cl)
  bb = Sys.time()
  print(paste0('KS Detection', 'c' = c, '---------------End, Time = ', bb - aa))
}


