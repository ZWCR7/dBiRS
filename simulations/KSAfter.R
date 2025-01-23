cv = c(0.25, 0.30)
alpha = 0.05

BIMList = read_bim('Gene_Sample/genotype.bim')

#chr = 1
BP = BIMList$pos
#RSID = BIMList$id

for (c in cv)
{
  RES_KS = c()
  for (s in 1:100)
  {
    begid = 1
    result.single = c(); result.window = c()
    
    n_snp = rep(0, 25)
    for (part in 1:25)
    {
      load(paste0('result_ks_b1/KS_c = ', c, 's = ', s, 'part', part, '.RData'))
      KSRES = result.fit
      
      result.single = rbind(result.single, KSRES$result.single)
      result.window = rbind(result.window, KSRES$result.window)
      
      n_snp[part] = nrow(KSRES$variant.info)
    }
    
    result.summary = as.data.frame(KS.summary(result.single, result.window, 5))
    indsel = which(result.summary$Qvalue < alpha)
    result.significant = result.summary[indsel, ]

    startind = rep(0, length(indsel)); endind = rep(0, length(indsel))
    if (nrow(result.significant) != 0)
    {
      for (l in 1:nrow(result.significant))
      {
        startind[l] = which(BP == result.significant$actual_start[l])
        endind[l] = which(BP == result.significant$actual_end[l])
      }

      res_KS = data.frame(startind, endind)
    }
    else
    {
      res_KS = 'There is no significant regions.'
    }
    
    RES_KS = c(RES_KS, res_KS = list(res_KS))
  }
  
  save(RES_KS, file = paste0('Gene_RES1/KS_COSI2_Lin_Signal_Binary_c = ', c, '.RData'))
}

