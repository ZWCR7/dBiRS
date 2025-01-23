library(genio)

nsimu = 100
num_true = 4

BIMList = read_bim('Gene_Sample/genotype.bim')
chr = 1
BP = BIMList$pos; index_trunc = max(which(BP < 5000*1000))
BP = BP[1:index_trunc]

p = length(BP)

EST_PROB_BBiRS = matrix(0, p, 2)
EST_PROB_KS = matrix(0, p, 2)
EST_PROB_QSCAN = matrix(0, p, 2)

RTrateBBiRS = matrix(0, nsimu, 2)
RFrateBBiRS1 = matrix(0, nsimu, 2)
RFrateBBiRS2 = matrix(0, nsimu, 2)
RFrateBBiRS3 = matrix(0, nsimu, 2)
PTrateBBiRS = matrix(0, nsimu, 2)
PFrateBBiRS1 = matrix(0, nsimu, 2)
PFrateBBiRS2 = matrix(0, nsimu, 2)
PFrateBBiRS3 = matrix(0, nsimu, 2)
RITPRBBiRS = matrix(0, num_true, 2)
ITPRBBiRS = matrix(0, num_true, 2)
IFDRBBiRS1 = matrix(0, num_true, 2)
IFDRBBiRS2 = matrix(0, num_true, 2)
IFDRBBiRS3 = matrix(0, num_true, 2)
LowerBBiRS = matrix(0, num_true, 2)
UpperBBiRS = matrix(0, num_true, 2)
PowerBBiRS = rep(nsimu, 2)

RTrateKS = matrix(0, nsimu, 2)
RFrateKS1 = matrix(0, nsimu, 2)
RFrateKS2 = matrix(0, nsimu, 2)
RFrateKS3 = matrix(0, nsimu, 2)
PTrateKS = matrix(0, nsimu, 2)
PFrateKS1 = matrix(0, nsimu, 2)
PFrateKS2 = matrix(0, nsimu, 2)
PFrateKS3 = matrix(0, nsimu, 2)
RITPRKS = matrix(0, num_true, 2)
ITPRKS = matrix(0, num_true, 2)
IFDRKS1 = matrix(0, num_true, 2)
IFDRKS2 = matrix(0, num_true, 2)
IFDRKS3 = matrix(0, num_true, 2)
LowerKS = matrix(0, num_true, 2)
UpperKS = matrix(0, num_true, 2)
PowerKS = rep(nsimu, 2)


RTrateQSCAN = matrix(0, nsimu, 2)
RFrateQSCAN1 = matrix(0, nsimu, 2)
RFrateQSCAN2 = matrix(0, nsimu, 2)
RFrateQSCAN3 = matrix(0, nsimu, 2)
PTrateQSCAN = matrix(0, nsimu, 2)
PFrateQSCAN1 = matrix(0, nsimu, 2)
PFrateQSCAN2 = matrix(0, nsimu, 2)
PFrateQSCAN3 = matrix(0, nsimu, 2)
RITPRQSCAN = matrix(0, num_true, 2)
ITPRQSCAN = matrix(0, num_true, 2)
IFDRQSCAN1 = matrix(0, num_true, 2)
IFDRQSCAN2 = matrix(0, num_true, 2)
IFDRQSCAN3 = matrix(0, num_true, 2)
LowerQSCAN = matrix(0, num_true, 2)
UpperQSCAN = matrix(0, num_true, 2)
PowerQSCAN = rep(nsimu, 2)

BIMList = read_bim('Gene_Sample/genotype.bim')
BP = BIMList$pos

#true.ind.frame = rbind(c(67513, 67876), c(218559, 218900), c(330769, 331125), c(455704, 456074), c(608297, 608650))
true.ind.frame = rbind(c(50165, 50504), c(124770, 125091), c(208781, 209123), c(291337, 291685))

Modify_True = function(d)
{
  true.ind.framedk = matrix(0, nrow(true.ind.frame), 2)
  true.inddk = c()
  
  for (i in 1:nrow(true.ind.frame))
  {
    true.start = BP[true.ind.frame[i, 1]]; true.end = BP[true.ind.frame[i, 2]]
    
    startdk = max(which(BP <= true.start - d*1000)) + 1; enddk = min(which(BP >= true.end + d*1000)) - 1
    true.ind.framedk[i, ] = c(startdk, enddk)
    true.inddk = c(true.inddk, startdk:enddk)
  }
  
  return(true.inddk)
}

pc = 0.1
covers_low = function(x) {return(quantile(x, pc, na.rm = T))}
covers_upper = function(x) {return(quantile(x, 1 - pc, na.rm = T))}

rdv = c(25, 50, 75); idv = c(1, 3, 5)

cv = c(0.12, 0.15)
for (setc in 1:length(cv))
{
  load(paste0('Gene_RES1/dBiRS_COSI2_Lin_Signal_Quantitive_c = ', cv[setc], '.RData'))
  load(paste0('Gene_RES1/QSCAN_COSI2_Lin_Signal_Quantitive_c = ', cv[setc], '.RData'))
  load(paste0('Gene_RES1/KS_COSI2_Lin_Signal_Quantitive_c = ', cv[setc], '.RData'))
  
  RTPRs = matrix(0, nsimu, num_true); PTPRs = matrix(0, nsimu, num_true); FDR1 = matrix(0, nsimu, num_true); FDR2 = matrix(0, nsimu, num_true); FDR3 = matrix(0, nsimu, num_true);
  Lowers = matrix(-1, nsimu, num_true); Uppers = matrix(-1, nsimu, num_true); EST_PROB = matrix(0, nsimu, p)
  
  for (i in 1:nsimu)
  {
    if (!is.character(RES_BBiRS[[i]]))
    {
      Resj = RES_BBiRS[[i]]$BSDCF_res
      num_detect = nrow(Resj)
      
      Rtrue.ind1 = Modify_True(rdv[1]); Rtrue.ind2 = Modify_True(rdv[2]); Rtrue.ind3 = Modify_True(rdv[3])
      Itrue.ind1 = Modify_True(idv[1]); Itrue.ind2 = Modify_True(idv[2]); Itrue.ind3 = Modify_True(idv[3])
      
      est.ind = c(); RFalse_num1 = 0;  RFalse_num2 = 0; RFalse_num3 = 0
      for (j in 1:num_detect)
      {
        est.indj = Resj[j, 1]:Resj[j, 2]
        est.ind = c(est.ind, est.indj)
        
        if (sum(est.indj %in% Rtrue.ind1) == 0) RFalse_num1 = RFalse_num1 + 1
        if (sum(est.indj %in% Rtrue.ind2) == 0) RFalse_num2 = RFalse_num2 + 1
        if (sum(est.indj %in% Rtrue.ind3) == 0) RFalse_num3 = RFalse_num3 + 1
      }
      
      est.ind = sort(unique(est.ind))
      EST_PROB[i, est.ind] = 1
      
      true.ind = c(); RTrue_num = 0
      for (k in 1:nrow(true.ind.frame))
      {
        true.indk = true.ind.frame[k, 1]:true.ind.frame[k, 2]
        true.ind = c(true.ind, true.indk)
        
        if (sum(true.indk %in% est.ind) > 0) RTrue_num = RTrue_num + 1
      }
      
      RTrateBBiRS[i, setc] = RTrue_num/nrow(true.ind.frame)
      RFrateBBiRS1[i, setc] = RFalse_num1/num_detect
      RFrateBBiRS2[i, setc] = RFalse_num2/num_detect
      RFrateBBiRS3[i, setc] = RFalse_num3/num_detect
      
      PTrue_num = sum(est.ind %in% true.ind)
      PFalse_num1 = length(est.ind) - sum(est.ind %in% Rtrue.ind1)
      PFalse_num2 = length(est.ind) - sum(est.ind %in% Rtrue.ind2)
      PFalse_num3 = length(est.ind) - sum(est.ind %in% Rtrue.ind3)
      PTrateBBiRS[i, setc] = PTrue_num/length(true.ind)
      PFrateBBiRS1[i, setc] = PFalse_num1/length(est.ind)
      PFrateBBiRS2[i, setc] = PFalse_num2/length(est.ind)
      PFrateBBiRS3[i, setc] = PFalse_num3/length(est.ind)
      
      for (k in 1:nrow(true.ind.frame))
      {
        true.indk = true.ind.frame[k, 1]:true.ind.frame[k, 2]
        est.indk = c()
        
        for (j in 1:num_detect)
        {
          est.indj = Resj[j, 1]:Resj[j, 2]
          if (sum(est.indj %in% true.indk) > 0)
          {
            est.indk = c(est.indk, est.indj)
          }
        }
        
        if (length(est.indk) > 0)
        {
          RTPRs[i, k] = 1
          
          k_True_num = sum(est.indk %in% true.indk)
          k_False_num1 = length(est.indk) - sum(est.indk %in% Itrue.ind1)
          k_False_num2 = length(est.indk) - sum(est.indk %in% Itrue.ind2)
          k_False_num3 = length(est.indk) - sum(est.indk %in% Itrue.ind3)
          
          PTPRs[i, k] = k_True_num/length(true.indk)
          FDR1[i, k] = k_False_num1/length(est.indk)
          FDR2[i, k] = k_False_num2/length(est.indk)
          FDR3[i, k] = k_False_num3/length(est.indk)
          
          Lowers[i, k] = min(est.indk); Uppers[i, k] = max(est.indk)
        }
      }
    }
    else
    {
      PowerBBiRS[setc] = PowerBBiRS[setc] - 1
    }
  }
  
  RITPRBBiRS[, setc] = apply(RTPRs, 2, mean)
  ITPRBBiRS[, setc] = apply(PTPRs, 2, mean)
  IFDRBBiRS1[, setc] = apply(FDR1, 2, mean)
  IFDRBBiRS2[, setc] = apply(FDR2, 2, mean)
  IFDRBBiRS3[, setc] = apply(FDR3, 2, mean)
  EST_PROB_BBiRS[, setc] = apply(EST_PROB, 2, mean)
  
  LowerBBiRS[, setc] = apply(Lowers, 2, covers_low)
  UpperBBiRS[, setc] = apply(Uppers, 2, covers_upper)
  ##################################################################################
  
  RTPRs = matrix(0, nsimu, num_true); PTPRs = matrix(0, nsimu, num_true); FDR1 = matrix(0, nsimu, num_true); FDR2 = matrix(0, nsimu, num_true); FDR3 = matrix(0, nsimu, num_true);
  Lowers = matrix(-1, nsimu, num_true); Uppers = matrix(-1, nsimu, num_true); EST_PROB = matrix(0, nsimu, p)
  
  for (i in 1:nsimu)
  {
    if (!is.character(RES_KS[[i]]))
    {
      Resj = RES_KS[[i]]
      num_detect = nrow(Resj)
      
      Rtrue.ind1 = Modify_True(rdv[1]); Rtrue.ind2 = Modify_True(rdv[2]); Rtrue.ind3 = Modify_True(rdv[3])
      Itrue.ind1 = Modify_True(idv[1]); Itrue.ind2 = Modify_True(idv[2]); Itrue.ind3 = Modify_True(idv[3])
      
      est.ind = c(); RFalse_num1 = 0;  RFalse_num2 = 0; RFalse_num3 = 0
      for (j in 1:num_detect)
      {
        est.indj = Resj[j, 1]:Resj[j, 2]
        est.ind = c(est.ind, est.indj)
        
        if (sum(est.indj %in% Rtrue.ind1) == 0) RFalse_num1 = RFalse_num1 + 1
        if (sum(est.indj %in% Rtrue.ind2) == 0) RFalse_num2 = RFalse_num2 + 1
        if (sum(est.indj %in% Rtrue.ind3) == 0) RFalse_num3 = RFalse_num3 + 1
      }
      
      est.ind = sort(unique(est.ind))
      EST_PROB[i, est.ind] = 1
      
      true.ind = c(); RTrue_num = 0
      for (k in 1:nrow(true.ind.frame))
      {
        true.indk = true.ind.frame[k, 1]:true.ind.frame[k, 2]
        true.ind = c(true.ind, true.indk)
        
        if (sum(true.indk %in% est.ind) > 0) RTrue_num = RTrue_num + 1
      }
      
      RTrateKS[i, setc] = RTrue_num/nrow(true.ind.frame)
      RFrateKS1[i, setc] = RFalse_num1/num_detect
      RFrateKS2[i, setc] = RFalse_num2/num_detect
      RFrateKS3[i, setc] = RFalse_num3/num_detect
      
      PTrue_num = sum(est.ind %in% true.ind)
      PFalse_num1 = length(est.ind) - sum(est.ind %in% Rtrue.ind1)
      PFalse_num2 = length(est.ind) - sum(est.ind %in% Rtrue.ind2)
      PFalse_num3 = length(est.ind) - sum(est.ind %in% Rtrue.ind3)
      PTrateKS[i, setc] = PTrue_num/length(true.ind)
      PFrateKS1[i, setc] = PFalse_num1/length(est.ind)
      PFrateKS2[i, setc] = PFalse_num2/length(est.ind)
      PFrateKS3[i, setc] = PFalse_num3/length(est.ind)
      
      for (k in 1:nrow(true.ind.frame))
      {
        true.indk = true.ind.frame[k, 1]:true.ind.frame[k, 2]
        est.indk = c()
        
        for (j in 1:num_detect)
        {
          est.indj = Resj[j, 1]:Resj[j, 2]
          if (sum(est.indj %in% true.indk) > 0)
          {
            est.indk = c(est.indk, est.indj)
          }
        }
        
        if (length(est.indk) > 0)
        {
          RTPRs[i, k] = 1
          
          k_True_num = sum(est.indk %in% true.indk)
          k_False_num1 = length(est.indk) - sum(est.indk %in% Itrue.ind1)
          k_False_num2 = length(est.indk) - sum(est.indk %in% Itrue.ind2)
          k_False_num3 = length(est.indk) - sum(est.indk %in% Itrue.ind3)
          
          PTPRs[i, k] = k_True_num/length(true.indk)
          FDR1[i, k] = k_False_num1/length(est.indk)
          FDR2[i, k] = k_False_num2/length(est.indk)
          FDR3[i, k] = k_False_num3/length(est.indk)
          
          Lowers[i, k] = min(est.indk); Uppers[i, k] = max(est.indk)
        }
      }
    }
    else
    {
      PowerKS[setc] = PowerKS[setc] - 1
    }
  }
  
  RITPRKS[, setc] = apply(RTPRs, 2, mean)
  ITPRKS[, setc] = apply(PTPRs, 2, mean)
  IFDRKS1[, setc] = apply(FDR1, 2, mean)
  IFDRKS2[, setc] = apply(FDR2, 2, mean)
  IFDRKS3[, setc] = apply(FDR3, 2, mean)
  
  LowerKS[, setc] = apply(Lowers, 2, covers_low)
  UpperKS[, setc] = apply(Uppers, 2, covers_upper)
  EST_PROB_KS[, setc] = apply(EST_PROB, 2, mean)
  ###################################################################################
  
  RTPRs = matrix(0, nsimu, num_true); PTPRs = matrix(0, nsimu, num_true); FDR1 = matrix(0, nsimu, num_true); FDR2 = matrix(0, nsimu, num_true); FDR3 = matrix(0, nsimu, num_true);
  Lowers = matrix(-1, nsimu, num_true); Uppers = matrix(-1, nsimu, num_true); EST_PROB = matrix(0, nsimu, p)

  for (i in 1:nsimu)
  {
    Resj = as.matrix(RES_QSCAN[[i]]$SCAN_res)

    if (ncol(Resj) == 1) Resj = t(Resj)
    if (sum(Resj) != 1)
    {
      num_detect = nrow(Resj)

      Rtrue.ind1 = Modify_True(rdv[1]); Rtrue.ind2 = Modify_True(rdv[2]); Rtrue.ind3 = Modify_True(rdv[3])
      Itrue.ind1 = Modify_True(idv[1]); Itrue.ind2 = Modify_True(idv[2]); Itrue.ind3 = Modify_True(idv[3])

      est.ind = c(); RFalse_num1 = 0;  RFalse_num2 = 0; RFalse_num3 = 0
      for (j in 1:num_detect)
      {
        est.indj = Resj[j, 2]:Resj[j, 3]
        est.ind = c(est.ind, est.indj)

        if (sum(est.indj %in% Rtrue.ind1) == 0) RFalse_num1 = RFalse_num1 + 1
        if (sum(est.indj %in% Rtrue.ind2) == 0) RFalse_num2 = RFalse_num2 + 1
        if (sum(est.indj %in% Rtrue.ind3) == 0) RFalse_num3 = RFalse_num3 + 1
      }

      est.ind = sort(unique(est.ind))
      EST_PROB[i, est.ind] = 1

      true.ind = c(); RTrue_num = 0
      for (k in 1:nrow(true.ind.frame))
      {
        true.indk = true.ind.frame[k, 1]:true.ind.frame[k, 2]
        true.ind = c(true.ind, true.indk)

        if (sum(true.indk %in% est.ind) > 0) RTrue_num = RTrue_num + 1
      }

      RTrateQSCAN[i, setc] = RTrue_num/nrow(true.ind.frame)
      RFrateQSCAN1[i, setc] = RFalse_num1/num_detect
      RFrateQSCAN2[i, setc] = RFalse_num2/num_detect
      RFrateQSCAN3[i, setc] = RFalse_num3/num_detect

      PTrue_num = sum(est.ind %in% true.ind)
      PFalse_num1 = length(est.ind) - sum(est.ind %in% Rtrue.ind1)
      PFalse_num2 = length(est.ind) - sum(est.ind %in% Rtrue.ind2)
      PFalse_num3 = length(est.ind) - sum(est.ind %in% Rtrue.ind3)
      PTrateQSCAN[i, setc] = PTrue_num/length(true.ind)
      PFrateQSCAN1[i, setc] = PFalse_num1/length(est.ind)
      PFrateQSCAN2[i, setc] = PFalse_num2/length(est.ind)
      PFrateQSCAN3[i, setc] = PFalse_num3/length(est.ind)

      for (k in 1:nrow(true.ind.frame))
      {
        true.indk = true.ind.frame[k, 1]:true.ind.frame[k, 2]
        est.indk = c()

        for (j in 1:num_detect)
        {
          est.indj = Resj[j, 2]:Resj[j, 3]
          if (sum(est.indj %in% true.indk) > 0)
          {
            est.indk = c(est.indk, est.indj)
          }
        }

        if (length(est.indk) > 0)
        {
          RTPRs[i, k] = 1

          k_True_num = sum(est.indk %in% true.indk)
          k_False_num1 = length(est.indk) - sum(est.indk %in% Itrue.ind1)
          k_False_num2 = length(est.indk) - sum(est.indk %in% Itrue.ind2)
          k_False_num3 = length(est.indk) - sum(est.indk %in% Itrue.ind3)

          PTPRs[i, k] = k_True_num/length(true.indk)
          FDR1[i, k] = k_False_num1/length(est.indk)
          FDR2[i, k] = k_False_num2/length(est.indk)
          FDR3[i, k] = k_False_num3/length(est.indk)

          Lowers[i, k] = min(est.indk); Uppers[i, k] = max(est.indk)
        }
      }
    }
    else
    {
      PowerQSCAN[setc] = PowerQSCAN[setc] - 1
    }
  }

  RITPRQSCAN[, setc] = apply(RTPRs, 2, mean)
  ITPRQSCAN[, setc] = apply(PTPRs, 2, mean)
  IFDRQSCAN1[, setc] = apply(FDR1, 2, mean)
  IFDRQSCAN2[, setc] = apply(FDR2, 2, mean)
  IFDRQSCAN3[, setc] = apply(FDR3, 2, mean)

  LowerQSCAN[, setc] = apply(Lowers, 2, covers_low)
  UpperQSCAN[, setc] = apply(Uppers, 2, covers_upper)
  EST_PROB_QSCAN[, setc] = apply(EST_PROB, 2, mean)
}

RTPRBBiRS = apply(RTrateBBiRS, 2, mean)
RFDRBBiRS1 = apply(RFrateBBiRS1, 2, mean)
RFDRBBiRS2 = apply(RFrateBBiRS2, 2, mean)
RFDRBBiRS3 = apply(RFrateBBiRS3, 2, mean)

PTPRBBiRS = apply(PTrateBBiRS, 2, mean)
PFDRBBiRS1 = apply(PFrateBBiRS1, 2, mean)
PFDRBBiRS2 = apply(PFrateBBiRS2, 2, mean)
PFDRBBiRS3 = apply(PFrateBBiRS3, 2, mean)

RTPRKS = apply(RTrateKS, 2, mean)
RFDRKS1 = apply(RFrateKS1, 2, mean)
RFDRKS2 = apply(RFrateKS2, 2, mean)
RFDRKS3 = apply(RFrateKS3, 2, mean)

PTPRKS = apply(PTrateKS, 2, mean)
PFDRKS1 = apply(PFrateKS1, 2, mean)
PFDRKS2 = apply(PFrateKS2, 2, mean)
PFDRKS3 = apply(PFrateKS3, 2, mean)

RTPRQSCAN = apply(RTrateQSCAN, 2, mean)
RFDRQSCAN1 = apply(RFrateQSCAN1, 2, mean)
RFDRQSCAN2 = apply(RFrateQSCAN2, 2, mean)
RFDRQSCAN3 = apply(RFrateQSCAN3, 2, mean)

PTPRQSCAN = apply(PTrateQSCAN, 2, mean)
PFDRQSCAN1 = apply(PFrateQSCAN1, 2, mean)
PFDRQSCAN2 = apply(PFrateQSCAN2, 2, mean)
PFDRQSCAN3 = apply(PFrateQSCAN3, 2, mean)

VRTPRBBiRS = apply(RTrateBBiRS, 2, sd)
VRFDRBBiRS1 = apply(RFrateBBiRS1, 2, sd)
VRFDRBBiRS2 = apply(RFrateBBiRS2, 2, sd)
VRFDRBBiRS3 = apply(RFrateBBiRS3, 2, sd)

VPTPRBBiRS = apply(PTrateBBiRS, 2, sd)
VPFDRBBiRS1 = apply(PFrateBBiRS1, 2, sd)
VPFDRBBiRS2 = apply(PFrateBBiRS2, 2, sd)
VPFDRBBiRS3 = apply(PFrateBBiRS3, 2, sd)

VRTPRKS = apply(RTrateKS, 2, sd)
VRFDRKS1 = apply(RFrateKS1, 2, sd)
VRFDRKS2 = apply(RFrateKS2, 2, sd)
VRFDRKS3 = apply(RFrateKS3, 2, sd)

VPTPRKS = apply(PTrateKS, 2, sd)
VPFDRKS1 = apply(PFrateKS1, 2, sd)
VPFDRKS2 = apply(PFrateKS2, 2, sd)
VPFDRKS3 = apply(PFrateKS3, 2, sd)

VRTPRQSCAN = apply(RTrateQSCAN, 2, sd)
VRFDRQSCAN1 = apply(RFrateQSCAN1, 2, sd)
VRFDRQSCAN2 = apply(RFrateQSCAN2, 2, sd)
VRFDRQSCAN3 = apply(RFrateQSCAN3, 2, sd)

VPTPRQSCAN = apply(PTrateQSCAN, 2, sd)
VPFDRQSCAN1 = apply(PFrateQSCAN1, 2, sd)
VPFDRQSCAN2 = apply(PFrateQSCAN2, 2, sd)
VPFDRQSCAN3 = apply(PFrateQSCAN3, 2, sd)

save.image(file = 'Gene_RES1/Summary_Simulation_Quantitive.RData')