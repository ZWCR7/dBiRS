CHR = 1

library(BiRS)
library(genio)
library(Matrix)
library(seqminer)
library(foreach)
library(doParallel)

source('distributed_BiRS.R')
load('AD20016_Sample.RData')
load(paste0('Blocks-50kb/Block-Chr', CHR, '.RData'))

partnum = nrow(Block_Index)

GenoPrefix = paste0('/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/ukb23158_c', CHR, '_b0_v1')
Null_GLM = NULL_GLM(FI_phenotype, FI_covariate, family = 'poisson')
#load('Null_GLM.RData')

Parallel_BBiRS = function(part)
{
  RES = distributed_BiRS(phenotype = FI_phenotype, X = FI_covariate, Null_GLM = Null_GLM, GenoPrefix = GenoPrefix, Block.ind = Block_Index[part, ], 
                   IfScale = F, subsample = sample_index, trunc = 3, alpha = 0.05)
  save(list = c('RES'), file = paste0('FI_results', CHR, '/part', part, '.RData'))
}

TT = ceiling(partnum/25)
for (tt in 1:TT)
{
  aa = Sys.time()
  start_tt = 25*(tt - 1) + 1
  end_tt = min(25*tt, partnum)
  
  coresNum = end_tt - start_tt + 1
  
  cl = makeCluster(coresNum)
  registerDoParallel(cl)
  
  RES = foreach(part = start_tt:end_tt, .packages = c('genio', 'BiRS', 'seqminer', 'Matrix')) %dopar% Parallel_BBiRS(part)
  
  stopImplicitCluster()
  stopCluster(cl)
  bb = Sys.time()
  
  print(paste0('Big part', tt, 'time = ', bb - aa))
}

system(paste0('tar -czvf FI_results', CHR, '.tar.gz FI_results', CHR))
system(paste0('dx upload FI_results', CHR, '.tar.gz'))


