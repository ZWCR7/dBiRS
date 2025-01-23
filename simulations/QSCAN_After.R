cv = c(0.12, 0.15)
alpha = 0.05
f = 0

for (c in cv)
{
  RES_QSCAN = c()
  for (s in 1:100)
  {
    begid = 1
    res = c(); resmost = c(); L20 = matrix(0, 1000, 100)
    
    n_snp = rep(0, 100)
    for (part in 1:100)
    {
      load(paste0('result_qscan_c1/QSCAN_c = ', c, 's = ', s, 'part', part, '.RData'))
      
      res_p = res_QSCAN$res$res
      resmost_p = res_QSCAN$resmost
      L20[, part] = res_QSCAN$L20
      
      if (sum(res) != 0)
      {
        res_p[, 2] = res_p[, 2] + begid - 1
        res_p[, 3] = res_p[, 3] + begid - 1
      }
      
      resmost_p[, 2] = resmost_p[, 2] + begid - 1
      resmost_p[, 3] = resmost_p[, 3] + begid - 1
      
      res = rbind(res, res_p)
      resmost = rbind(resmost, resmost_p)
      
      begid = begid + res_QSCAN$n_snp
      
      n_snp[part] = res_QSCAN$n_snp
    }
    
    ThresV = apply(L20, 1, max)
    ThresMax = quantile(ThresV, 1 - alpha)
    
    res = res[res[, 1] > ThresMax, ]
    
    if(length(res) == 0)
    {
      res = c(0, 0, 0, 1)
    }
    
    if(length(res) > 4)
    {
      res = regionfilter(res, f)
      if(length(res) == 4)
      {
        res[4] = mean(ThresV > res[1])
      }
      else
      {
        res[,4] = apply(res, 1, function(z) mean(ThresV > z[1]))
      }
    }
    
    mostnum = which.max(resmost[, 1])
    resmost = resmost[mostnum, ]
    resmost[4] = mean(ThresV > resmost[1])
    
    res_qscan = list(SCAN_res = res, SCAN_top1 = resmost, SCAN_thres = ThresMax, n_snp = n_snp)
    
    RES_QSCAN = c(RES_QSCAN, res_QSCAN = list(res_qscan))
  }
  
  save(RES_QSCAN, file = paste0('Gene_RES1/QSCAN_COSI2_Lin_Signal_Quantitive_c = ', c, '.RData'))
}

