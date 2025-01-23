library(BiRS)

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
