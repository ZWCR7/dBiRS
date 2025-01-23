BiRS_Generator = function(phenotype, X, G, mu0, working, sigma, fam, foldlen, trunc, MB = 1000, alpha = 0.05, ReMax = 10, bigmemory = T)
{
  n = length(phenotype); p = ncol(G)
  
  if (bigmemory)
  {
    SumStat = SummaryGenerator(G, phenotype, X, mu0, working, sigma, fam, MB)
  }
  else
  {
    SumStat = SummaryGenerator_M(G, phenotype, X, mu0, working, sigma, fam, MB)
  }
  
  UnVec = abs(SumStat$UnV)/sqrt(n); SneX = abs(SumStat$SneX)/sqrt(n)
  rm(SumStat); gc()
  
  UnStat = UnVec
  ThresStat = SneX
  
  quantileMax = apply(SneX, 1, max)
  bmax_prm = quantile(quantileMax, 1 - alpha)
  
  Signal = rep(0, p)
  split = ceiling(p/foldlen)
  
  S.split = list()
  
  if (split == 1)
  {
    S.split = c(S.split, list(1:p))
  }
  if (split != 1)
  {
    for (i in 1:(split - 1))
    {
      S.split = c(S.split, list(((i-1)*foldlen+1):(i*foldlen)))
    }
    S.split = c(S.split, list(((split-1)*foldlen+1):p))
  }
  
  Unprm = rep(0, split)
  for (sp in 1:split)
  {
    Unprm[sp] = max(UnVec[S.split[[sp]]])
  }
  
  reject_ind = rep(0, split)
  Slist = list()
  for (sp in 1:split)
  {
    if (Unprm[sp] > bmax_prm)
    {
      reject_ind[sp] = 1
      ps = length(S.split[[sp]]); ms = ceiling(ps/2)
      Slist = c(Slist, list(S.split[[sp]][1:ms]))
      Slist = c(Slist, list(S.split[[sp]][(ms + 1):ps]))
    }
  }
  
  if (sum(reject_ind) == 0) return(list(BSDCF = NULL, quantileMax = quantileMax, quantileMin = rep(-9999, MB)))
  
  Signal = BiSearchDCF(UnVec, SneX, Slist, Signal, trunc, MB, alpha)
  
  ind = which(Signal == 1)
  UnVec[ind] = 0
  SneX[ ,ind] = 0
  loop.res = Rejec(UnVec, SneX, alpha)
  loop.ind = loop.res$rej
  loop.bound = loop.res$bd
  
  rm(loop.res)
  
  index = 0
  while ((loop.ind == T) & (index <= ReMax))
  {
    Unprm.loop = rep(0, split)
    for (sp in 1:split)
    {
      Unprm.loop[sp] = max(UnVec[S.split[[sp]]])
    }
    
    Slist.loop = list()
    for (sp in 1:split)
    {
      if (Unprm.loop[sp] > loop.bound)
      {
        p.loop = length(S.split[[sp]]); m.loop = ceiling(p.loop/2)
        Slist.loop = c(Slist.loop, list(S.split[[sp]][1:m.loop]))
        Slist.loop = c(Slist.loop, list(S.split[[sp]][(m.loop + 1):p.loop]))
      }
    }
    
    Signal = BiSearchDCF(UnVec, SneX, Slist.loop, Signal, trunc, MB, alpha)
    
    ind = which(Signal == 1)
    UnVec[ind] = 0
    SneX[ ,ind] = 0
    loop.res = Rejec(UnVec, SneX, alpha)
    loop.ind = loop.res$rej
    loop.bound = loop.res$bd
    
    rm(loop.res)
    
    index = index + 1
  }
  
  indF = which(Signal == 1)
  #UnSelect = UnStat[indF]
  quantileMin = apply(ThresStat[, indF], 1, max)
  
  nF = length(indF)
  startind = indF[1]; endind = NULL
  for (i in 1:(nF - 1))
  {
    if ((indF[i + 1] - indF[i]) > 1)
    {
      endind = c(endind, indF[i])
      startind = c(startind, indF[i + 1])
    }
  }
  endind = c(endind, indF[nF])
  
  Tn = length(startind)
  for (i in 1:length(startind))
  {
    Tn[i] = max(UnStat[startind[i]:endind[i]])
  }
  #MaxStat = rep(UnMax, nF)
  
  BSDCF = data.frame(startind, endind, Tn)
  
  return(list(BSDCF = BSDCF, quantileMax = quantileMax, quantileMin = quantileMin))
}

Impute = function(Z, impute.method)
{
  p = dim(Z)[2]
  
  if(impute.method == "random")
  {
    for(i in 1:p)
    {
      IDX = which(is.na(Z[,i]))
      if(length(IDX) > 0)
      {
        maf1 = mean(Z[-IDX,i])/2
        Z[IDX,i] = rbinom(length(IDX),2,maf1)
      }
    }
  }
  
  if(impute.method == "fixed")
  {
    for(i in 1:p)
    {
      IDX = which(is.na(Z[,i]))
      if(length(IDX) > 0)
      {
        maf1 = mean(Z[-IDX,i])/2
        Z[IDX,i] = 2*maf1
      }
    }
  }
  
  if(impute.method == "bestguess") 
  {
    for(i in 1:p)
    {
      IDX = which(is.na(Z[,i]))
      if(length(IDX) > 0)
      {
        maf1 = mean(Z[-IDX,i])/2
        Z[IDX,i] = round(2*maf1)
      }
    }
  }  
  return(as.matrix(Z))
}

distributed_BiRS = function(phenotype, X, Null_GLM, GenoPrefix, Block.ind, IfScale = T, subsample, trunc, alpha = 0.05, MB = 1000, ReMax = 10, bigmemory = T)
{
  n = length(phenotype)
  
  mu0 = Null_GLM$mu0; working = Null_GLM$working; sigma = Null_GLM$sigma; fam = Null_GLM$fam
  
  chrb = Block.ind[1]
  startb = Block.ind[2]; endb = Block.ind[3]
  Gb = -1*(readPlinkToMatrixByIndex(plinkFilePrefix = GenoPrefix, sampleIndex = subsample, markerIndex = startb:endb) - 2)
  
  Gb[Gb < 0 | Gb > 2] = NA
  
  MISS.freq = colMeans(is.na(Gb))
  index.miss = which(MISS.freq > 0.05)
  
  if (length(index.miss) > 0)
  {
    Gb = Gb[, -index.miss]
  }
  
  if (sum(is.na(Gb)) > 0)
  {
    Gb = Impute(Gb, 'fixed')
  }
  
  MACs = colSums(Gb)
  index.ultra = which(MACs <= 1)
  
  if (length(index.ultra) > 0)
  {
    Gb = Gb[, -index.ultra]
  }
  
  #Gb = scale(Gb, center = F, scale = T)
  foldlen = ncol(Gb); nsnp = ncol(Gb); mafs = colMeans(Gb)/2
  Gb = scale(Gb, center = T, scale = IfScale)
  
  BiRSRES = BiRS_Generator(phenotype, X, Gb, mu0, working, sigma, fam, foldlen, trunc, MB, alpha, ReMax, bigmemory)
  
  return(list(BiRSRES = BiRSRES, nsnp = nsnp, chr = chrb, bp = colnames(Gb), mafs = mafs))
}

NULL_GLM = function(phenotype, X, family)
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
  
  return(list(sigma = sigma, fam = fam, working = working, mu0 = mu0))
}