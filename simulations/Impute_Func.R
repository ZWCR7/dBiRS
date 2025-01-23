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

Impute_v = function(Z, impute.method)
{
  p = length(Z)
  
  if(impute.method == "random")
  {
    IDX = which(is.na(Z))
    if(length(IDX) > 0)
    {
      maf1 = mean(Z[-IDX])/2
      Z[IDX] = rbinom(length(IDX),2,maf1)
    }
  }
  
  if(impute.method == "fixed")
  {
    IDX = which(is.na(Z))
    if(length(IDX) > 0)
    {
      maf1 = mean(Z[-IDX])/2
      Z[IDX] = 2*maf1
    }
  }
  
  if(impute.method == "bestguess") 
  {
    IDX = which(is.na(Z))
    if(length(IDX) > 0)
    {
      maf1 = mean(Z[-IDX])/2
      Z[IDX] = round(2*maf1)
    }
  }  
  return(Z)
}
