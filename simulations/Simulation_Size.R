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
source('DistributedBiRS.R')
source('UniBlockQSCAN.R')

inv_logi = function(x) exp(x)/(1 + exp(x))

alphav = c(0.05)
for (alpha in alphav)
{
  Ycore = X %*% gamma

  nsimu = 1000

  Simu_dBiRS = function(s)
  {
    set.seed(s)

    epsilon = rnorm(n, 0, 1)
    Y = Ycore + epsilon
    phenotype = as.vector(Y)

    XTrans = scale(X, center = T, scale = F)

    QS_prefix = QSCAN_Prefix(phenotype = phenotype, X = X, family = 'gaussian')
    working = QS_prefix$working; sigma = QS_prefix$sigma; fam = QS_prefix$fam; mu0 = QS_prefix$mu0

    for (part in 1:nrow(Block.ind))
    {
      ind.part = Block.ind[part, ]
      res_dBiRS = dBiRS_Generator(phenotype = phenotype, X = XTrans, GenoPrefix = GenoPrefix, Block.ind = ind.part,
                                  sample_index = select_sample, IfScale = F, family = 'gaussian', trunc = 6, MB = 1000, alpha = alpha)
      
      res_QSCAN = QSCAN_Block(phenotype = phenotype, X = X, GenoPrefix = GenoPrefix, Block = ind.part, sample_index = select_sample,
                              working = working, sigma = sigma, fam = fam, mu0 = mu0, Lmax = 350, Lmin = 320)

      save(list = c('res_dBiRS'), file = paste0('result_dbirs_c_size/dBiRS_alpha = ', alpha, 's = ', s, 'part', part, '.RData'))
      save(list = c('res_QSCAN'), file = paste0('result_qscan_c_size/QSCAN_alpha = ', alpha, 's = ', s, 'part', part, '.RData'))
    }
  }

  print('Block BiRS Detection-----------------Start')
  aa = Sys.time()
  cl = makeCluster(100)
  registerDoParallel(cl)

  RES_BBiRS = foreach(s = 1:nsimu, .packages = c("BiRS", "MASS", "QSCAN", "Matrix", "seqminer", "genio", "KnockoffScreen")) %dopar% Simu_dBiRS(s)

  stopImplicitCluster()
  stopCluster(cl)
  bb = Sys.time()
  print(paste0('Block BiRS Detection', 'alpha' = alpha, '---------------End, Time = ', bb - aa))
}

#############################################################################################################################################################################################################
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
source('DistributedBiRS.R')
source('UniBlockQSCAN.R')

inv_logi = function(x) exp(x)/(1 + exp(x))

alphav = c(0.05)
for (alpha in alphav)
{
  Ycore = X %*% gamma
  
  nsimu = 1000
  
  Simu_dBiRS = function(s)
  {
    set.seed(s)
    
    eta = Ycore
    pi = inv_logi(eta)
    Y = rbinom(n, size = 1, prob = pi)
    phenotype = as.vector(Y)
    
    XTrans = scale(X, center = T, scale = F)
    
    QS_prefix = QSCAN_Prefix(phenotype = phenotype, X = X, family = 'binomial')
    working = QS_prefix$working; sigma = QS_prefix$sigma; fam = QS_prefix$fam; mu0 = QS_prefix$mu0
    
    for (part in 1:nrow(Block.ind))
    {
      ind.part = Block.ind[part, ]
      res_dBiRS = dBiRS_Generator(phenotype = phenotype, X = XTrans, GenoPrefix = GenoPrefix, Block.ind = ind.part,
                                  sample_index = select_sample, IfScale = F, family = 'binomial', trunc = 6, MB = 1000, alpha = alpha)
      
      res_QSCAN = QSCAN_Block(phenotype = phenotype, X = X, GenoPrefix = GenoPrefix, Block = ind.part, sample_index = select_sample,
                              working = working, sigma = sigma, fam = fam, mu0 = mu0, Lmax = 350, Lmin = 320)
      
      save(list = c('res_dBiRS'), file = paste0('result_dbirs_b_size/dBiRS_alpha = ', alpha, 's = ', s, 'part', part, '.RData'))
      save(list = c('res_QSCAN'), file = paste0('result_qscan_b_size/QSCAN_alpha = ', alpha, 's = ', s, 'part', part, '.RData'))
    }
  }
  
  print('Block BiRS Detection-----------------Start')
  aa = Sys.time()
  cl = makeCluster(100)
  registerDoParallel(cl)
  
  RES_BBiRS = foreach(s = 1:nsimu, .packages = c("BiRS", "MASS", "QSCAN", "Matrix", "seqminer", "genio", "KnockoffScreen")) %dopar% Simu_dBiRS(s)
  
  stopImplicitCluster()
  stopCluster(cl)
  bb = Sys.time()
  print(paste0('Block BiRS Detection', 'alpha' = alpha, '---------------End, Time = ', bb - aa))
}

