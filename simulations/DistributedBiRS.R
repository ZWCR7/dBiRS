library(BiRS)
library(genio)
library(Matrix)
library(seqminer)

source('BiRS-Generator.R')
source('Impute_Func.R')

dBiRS_Generator = function(phenotype, X, GenoPrefix, Block.ind, sample_index, IfScale = T, family, trunc, MB = 1000, alpha = 0.05, ReMax = 10, bigmemory = T)
{
  n = length(phenotype)
  if (family == "gaussian")
  {
    lmnull = lm(phenotype ~ X)
    sigma = summary(lmnull)$sigma
    fam = 0
    working = rep(1, n)
    mu0 = lmnull$fitted
  }
  
  if (family != "gaussian")
  {
    glmnull = glm(phenotype ~ X, family = family)
    sigma = sqrt(summary(glmnull)$dispersion)
    fam = 1
    working = glmnull$weights
    mu0 = glmnull$fitted
  }
  
  chrb = Block.ind[1]; startb = Block.ind[2]; endb = Block.ind[3]
  
  Gb = -1*(readPlinkToMatrixByIndex(plinkFilePrefix = GenoPrefix, sampleIndex = sample_index, markerIndex = startb:endb) - 2)
  
  Gb[Gb < 0 | Gb > 2] = NA
  if (sum(is.na(Gb)) > 0)
  {
    Gb = apply(Gb, 2, Impute_v, impute.method = 'bestguess')
  }
  
  
  foldlen = ncol(Gb); n_snp = ncol(Gb)
  Gb = scale(Gb, center = T, scale = IfScale)
  
  BiRSRES = BiRS_Generator(phenotype, X, Gb, mu0, working, sigma, fam, foldlen, trunc, MB, alpha, ReMax, bigmemory)
  
  return(list(BiRSRES = BiRSRES, n_snp = n_snp))
}