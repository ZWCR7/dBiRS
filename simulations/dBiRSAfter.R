cv = c(0.12, 0.15)
alpha = 0.05

SigDetDCFGe = function(UnVec, SneX, foldlen = length(UnVec), trunc = 0, alpha, ReMax = 10)
{
  p = length(UnVec)
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
  
  SneX.max = apply(SneX, 1, max)
  bmax_prm = quantile(SneX.max, prob = 1 - alpha)
  
  UFC = UnVec
  
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
  
  if (sum(reject_ind) == 0) return("there is no evidence that there exist signal region")
  
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
  
  
  return(Signal)
}

for (c in cv)
{
  RES_BBiRS = c()
  
  for (s in 1:100)
  {
    UnVec = NULL; SneX = NULL; partinf = NULL
    start_index = 0
    
    n_snp = rep(0, 100)
    for (part in 1:100)
    {
      load(paste0('result_dbirs_c1/dBiRS_c = ', c, 's = ', s, 'part', part, '.RData'))
      
      SneX = cbind(SneX, res_dBiRS$BiRSRES$quantileMax)
      #quantileMin = pmax(quantileMin, RES$BiRSRES$quantileMin)
      
      if (!is.null(res_dBiRS$BiRSRES$BSDCF))
      {
        UnVec = c(UnVec, max(abs(res_dBiRS$BiRSRES$BSDCF$Tn)))
      }
      else
      {
        UnVec = c(UnVec, 0)
      }
    }
    
    Block_Signal = SigDetDCFGe(UnVec = UnVec, SneX = SneX, alpha = alpha)
    significant_Block = which(Block_Signal == 1)
    
    if (length(significant_Block > 0))
    {
      quantileMin = rep(-9999, 1000)
      
      BSDCF_res = NULL; block_index = NULL
      for (part in 1:100)
      {
        load(paste0('result_dbirs_c1/dBiRS_c = ', c, 's = ', s, 'part', part, '.RData'))
        
        if (sum(significant_Block == part) > 0)
        {
          quantileMin = pmax(quantileMin, res_dBiRS$BiRSRES$quantileMin)
          BSDCFb = res_dBiRS$BiRSRES$BSDCF
          
          startind = NULL; endind = NULL 
          
          BSDCFb[, 1] = BSDCFb[, 1] + start_index
          BSDCFb[, 2] = BSDCFb[, 2] + start_index
          
          BSDCFb = cbind(BSDCFb, rep(part, nrow(BSDCFb)))
          
          BSDCF_res = rbind(BSDCF_res, BSDCFb)
          block_index = c(block_index, part)
        }
        
        start_index = start_index + res_dBiRS$n_snp
      }
      
      thresMin = quantile(quantileMin, 1 - alpha)
      
      index_thres = which(BSDCF_res[, 3] > thresMin)
      BSDCF_res = BSDCF_res[index_thres, ]
      
      res_BBiRS = list(BSDCF_res = BSDCF_res)
    }
    else
    {
      res_BBiRS = 'There is no significant regions.'
    }
    
    RES_BBiRS = c(RES_BBiRS, res_BBiRS = list(res_BBiRS))
  }
  
  save(RES_BBiRS, file = paste0('Gene_RES1/dBiRS_COSI2_Lin_Signal_Quantitive_c = ', c, '.RData'))
}
