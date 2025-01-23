QSCAN_Block = function(phenotype, X, GenoPrefix, Block, sample_index, working, sigma, fam, mu0, Lmax, Lmin, steplength = 1, times = 1000, alpha = 0.05, f = 0)
{
  # B_NUM = nrow(Block.ind)
  # n_snp = rep(0, B_NUM)
  n = length(phenotype)
  begid = 1
  
  chrb = Block[1]
  startb = Block[2]; endb = Block[3]
  
  Gb = -1*(readPlinkToMatrixByIndex(plinkFilePrefix = GenoPrefix, sampleIndex = sample_index, markerIndex = startb:endb) - 2)
  
  Gb[Gb < 0 | Gb > 2] = NA
  if (sum(is.na(Gb)) > 0)
  {
    Gb = apply(Gb, 2, Impute_v, impute.method = 'bestguess')
  }
  
  maf = colMeans(Gb)/2
  
  #Gb = G[, startb:endb]
  #n_snp[b] = ncol(Gb)
  Gb = scale(Gb, center = T, scale = F)
  Gb = as(Gb, 'sparseMatrix')
  
  weights = dbeta(maf, 1, 1)
  emL20 = Q_SCAN_Thres(Gb, X, working, sigma, fam, times, Lmax, Lmin, weights)
  
  print('Thres---end')
  
  th0 = quantile(emL20, 1 - alpha)
  
  res = Q_SCAN_Search(Gb, X, working, sigma, fam, phenotype, mu0, th0, Lmax, Lmin, begid, f, weights)
  
  print('Search---end')
  resmost = res$resmost
  
  #res = rbind(res, restemp$res)
  #resmost = rbind(resmost, restemp$resmost)
  
  L20 = emL20
  
  n_snp = ncol(Gb)
  
  return(list(res = res, resmost = resmost, n_snp = n_snp, L20 = L20))
}

QSCAN_Prefix = function(phenotype, X, family)
{
  samplesize = length(phenotype)
  
  if(family == "gaussian")
  {
    lmnull = lm(phenotype ~ -1 + X)
    sigma = summary(lmnull)$sigma
    
    fam = 0
    working = rep(1, samplesize)
    
    mu0 = lmnull$fitted
  }
  
  if(family != "gaussian")
  {
    # fit global null model
    glmnull = glm(phenotype ~ -1 + X, family = family)
    sigma = sqrt(summary(glmnull)$dispersion)
    
    fam = 1
    working = glmnull$weights
    
    mu0 = glmnull$fitted
  }
  
  return(list(working = working, sigma = sigma, fam = fam, mu0 = mu0))
}