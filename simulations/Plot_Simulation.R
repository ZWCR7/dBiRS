rm(list = ls())
load('Gene_RES1/Summary_Simulation_Binary.RData')

library(ggplot2)
library(ggbreak)

load('true_signal_binary.RData')

BIMList = read_bim('Gene_Sample/genotype.bim')

chr = 1
BP = BIMList$pos; index_trunc = max(which(BP < 5000*1000))

BP = BP[1:index_trunc]
RSID = BIMList$id[1:index_trunc]

posmin = min(BP); posmax = max(BP)

pos.tag = seq(posmin, posmax, by = 5000/2)
window.bed = cbind(chr, pos.tag, pos.tag + 5000 - 3)
window.bed = window.bed[order(as.numeric(window.bed[, 2])),]

ind.start = 1; Block.bed = NULL; Block.ind = NULL
pos.end = 0

while(pos.end < posmax)
{
  pos.start = BP[ind.start]
  pos.trunc = pos.start + 50*1000
  
  ind.end = max(which(BP <= pos.trunc))
  pos.end = BP[ind.end]
  
  Block.bed = rbind(Block.bed, c(chr, pos.start, pos.end))
  Block.ind = rbind(Block.ind, c(chr, ind.start, ind.end))
  
  print(ind.end - ind.start)
  ind.start = ind.end + 1
}

locuind = c(15, 36, 60, 84)

for (setc in 1:length(cv))
{
  signals = abs(mu[, setc]) 
  
  cover_plot_dBiRS = ggplot()
  cover_plot_KS = ggplot()
  cover_plot_QSCAN = ggplot()
  
  X_dBiRS = which(EST_PROB_BBiRS[, setc] > 0.95)
  X_KS = which(EST_PROB_KS[, setc] > 0.95)
  X_QSCAN = which(EST_PROB_QSCAN[, setc] > 0.95)
  
  Y_dBiRS = signals[which(EST_PROB_BBiRS[, setc] > 0.95)]
  Y_KS = signals[which(EST_PROB_KS[, setc] > 0.95)]
  Y_QSCAN = signals[which(EST_PROB_QSCAN[, setc] > 0.95)]
  
  Data_dBiRS = data.frame(X_dBiRS, Y_dBiRS)
  Data_KS = data.frame(X_KS, Y_KS)
  Data_QSCAN = data.frame(X_QSCAN, Y_QSCAN)
  
  X_prob = 1:p
  prob_dBiRS = EST_PROB_BBiRS[, setc]
  Data_dBiRS_Prob = data.frame(X_prob, prob_dBiRS)
  prob_plot_dBiRS = ggplot() + geom_point(data = Data_dBiRS_Prob, aes(x = X_prob, y = prob_dBiRS), size = 2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Probability")
  
  png(filename = paste0('Figures/dBiRS_Binary_Probability_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(prob_plot_dBiRS)
  dev.off()
  
  X_prob = 1:p
  prob_KS = EST_PROB_KS[, setc]
  Data_KS_Prob = data.frame(X_prob, prob_KS)
  prob_plot_KS = ggplot() + geom_point(data = Data_KS_Prob, aes(x = X_prob, y = prob_KS), size = 2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") + 
    theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Probability")
  
  png(filename = paste0('Figures/KS_Binary_Probability_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(prob_plot_KS)
  dev.off()
  
  X_prob = 1:p
  prob_QSCAN = EST_PROB_QSCAN[, setc]
  Data_QSCAN_Prob = data.frame(X_prob, prob_QSCAN)
  prob_plot_QSCAN = ggplot() + geom_point(data = Data_QSCAN_Prob, aes(x = X_prob, y = prob_QSCAN), size = 2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Probability")
  
  png(filename = paste0('Figures/QSCAN_Binary_Probability_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(prob_plot_QSCAN)
  dev.off()
  
  
  breaks = c()
  for (i in 1:num_true)
  {
    start = Block.ind[locuind[i], 2]; end = Block.ind[locuind[i], 3]
    
    X = start:end; Y = signals[start:end]
    Data = data.frame(X, Y)
    
    cover_plot_dBiRS = cover_plot_dBiRS + geom_point(data = Data, aes(x = X, y = Y), size = 2, col = 'grey') + 
      annotate("rect", xmin = true.ind.frame[i, 1], xmax = true.ind.frame[i, 2], ymin = 0, ymax = Inf, fill= "pink", alpha = 0.4) 
    
    cover_plot_KS = cover_plot_KS + geom_point(data = Data, aes(x = X, y = Y), size = 2, col = 'grey') + 
      annotate("rect", xmin = true.ind.frame[i, 1], xmax = true.ind.frame[i, 2], ymin = 0, ymax = Inf, fill= "pink", alpha = 0.4) 
    
    cover_plot_QSCAN = cover_plot_QSCAN + geom_point(data = Data, aes(x = X, y = Y), size = 2, col = 'grey') + 
      annotate("rect", xmin = true.ind.frame[i, 1], xmax = true.ind.frame[i, 2], ymin = 0, ymax = Inf, fill= "pink", alpha = 0.4) 
    
    breaks = c(breaks, floor((true.ind.frame[i, 1] + true.ind.frame[i, 2])/2))
  }
  
  cover_plot_dBiRS = cover_plot_dBiRS + geom_point(data = Data_dBiRS, aes(x = X_dBiRS, y = Y_dBiRS), size = 2, col = 'red')
  cover_plot_KS = cover_plot_KS + geom_point(data = Data_KS, aes(x = X_KS, y = Y_KS), size = 2, col = 'blue')
  cover_plot_QSCAN = cover_plot_QSCAN + geom_point(data = Data_QSCAN, aes(x = X_QSCAN, y = Y_QSCAN), size = 2, col = 'green')
  
  cover_plot_dBiRS = cover_plot_dBiRS + theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) + 
    scale_y_continuous(breaks = NULL) +
    scale_x_break(c(Block.ind[locuind[1], 3], Block.ind[locuind[2], 2])) +
    scale_x_break(c(Block.ind[locuind[2], 3], Block.ind[locuind[3], 2])) +
    scale_x_break(c(Block.ind[locuind[3], 3], Block.ind[locuind[4], 2])) +
    scale_x_continuous(breaks = breaks) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Underlying Effect")
  
  cover_plot_QSCAN = cover_plot_QSCAN + theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) + 
    scale_y_continuous(breaks = NULL) +
    scale_x_break(c(Block.ind[locuind[1], 3], Block.ind[locuind[2], 2])) +
    scale_x_break(c(Block.ind[locuind[2], 3], Block.ind[locuind[3], 2])) +
    scale_x_break(c(Block.ind[locuind[3], 3], Block.ind[locuind[4], 2])) +
    scale_x_continuous(breaks = breaks) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Underlying Effect")
  
  cover_plot_KS = cover_plot_KS + theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) + 
    scale_y_continuous(breaks = NULL) +
    scale_x_break(c(Block.ind[locuind[1], 3], Block.ind[locuind[2], 2])) +
    scale_x_break(c(Block.ind[locuind[2], 3], Block.ind[locuind[3], 2])) +
    scale_x_break(c(Block.ind[locuind[3], 3], Block.ind[locuind[4], 2])) +
    scale_x_continuous(breaks = breaks) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Underlying Effect")
  
  png(filename = paste0('Figures/dBiRS_Binary_Coverage_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(cover_plot_dBiRS)
  dev.off()
  
  png(filename = paste0('Figures/QSCAN_Binary_Coverage_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(cover_plot_QSCAN)
  dev.off()
  
  png(filename = paste0('Figures/KS_Binary_Coverage_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(cover_plot_KS)
  dev.off()
}
##########################################################################################################################################################################

rm(list = ls())
load('Gene_RES1/Summary_Simulation_Quantitive.RData')

library(ggplot2)
library(ggbreak)

load('true_signal_quantitive.RData')

BIMList = read_bim('Gene_Sample/genotype.bim')

chr = 1
BP = BIMList$pos; index_trunc = max(which(BP < 5000*1000))

BP = BP[1:index_trunc]
RSID = BIMList$id[1:index_trunc]

posmin = min(BP); posmax = max(BP)

pos.tag = seq(posmin, posmax, by = 5000/2)
window.bed = cbind(chr, pos.tag, pos.tag + 5000 - 3)
window.bed = window.bed[order(as.numeric(window.bed[, 2])),]

ind.start = 1; Block.bed = NULL; Block.ind = NULL
pos.end = 0

while(pos.end < posmax)
{
  pos.start = BP[ind.start]
  pos.trunc = pos.start + 50*1000
  
  ind.end = max(which(BP <= pos.trunc))
  pos.end = BP[ind.end]
  
  Block.bed = rbind(Block.bed, c(chr, pos.start, pos.end))
  Block.ind = rbind(Block.ind, c(chr, ind.start, ind.end))
  
  print(ind.end - ind.start)
  ind.start = ind.end + 1
}

locuind = c(15, 36, 60, 84)

for (setc in 1:length(cv))
{
  signals = abs(mu[, setc]) 
  
  cover_plot_dBiRS = ggplot()
  cover_plot_KS = ggplot()
  cover_plot_QSCAN = ggplot()
  
  X_dBiRS = which(EST_PROB_BBiRS[, setc] > 0.95)
  X_KS = which(EST_PROB_KS[, setc] > 0.95)
  X_QSCAN = which(EST_PROB_QSCAN[, setc] > 0.95)
  
  Y_dBiRS = signals[which(EST_PROB_BBiRS[, setc] > 0.95)]
  Y_KS = signals[which(EST_PROB_KS[, setc] > 0.95)]
  Y_QSCAN = signals[which(EST_PROB_QSCAN[, setc] > 0.95)]
  
  Data_dBiRS = data.frame(X_dBiRS, Y_dBiRS)
  Data_KS = data.frame(X_KS, Y_KS)
  Data_QSCAN = data.frame(X_QSCAN, Y_QSCAN)
  
  X_prob = 1:p
  prob_dBiRS = EST_PROB_BBiRS[, setc]
  Data_dBiRS_Prob = data.frame(X_prob, prob_dBiRS)
  prob_plot_dBiRS = ggplot() + geom_point(data = Data_dBiRS_Prob, aes(x = X_prob, y = prob_dBiRS), size = 2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Probability")
  
  png(filename = paste0('Figures/dBiRS_Quantitive_Probability_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(prob_plot_dBiRS)
  dev.off()
  
  X_prob = 1:p
  prob_KS = EST_PROB_KS[, setc]
  Data_KS_Prob = data.frame(X_prob, prob_KS)
  prob_plot_KS = ggplot() + geom_point(data = Data_KS_Prob, aes(x = X_prob, y = prob_KS), size = 2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Probability")
  
  png(filename = paste0('Figures/KS_Quantitive_Probability_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(prob_plot_KS)
  dev.off()
  
  X_prob = 1:p
  prob_QSCAN = EST_PROB_QSCAN[, setc]
  Data_QSCAN_Prob = data.frame(X_prob, prob_QSCAN)
  prob_plot_QSCAN = ggplot() + geom_point(data = Data_QSCAN_Prob, aes(x = X_prob, y = prob_QSCAN), size = 2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Probability")
  
  png(filename = paste0('Figures/QSCAN_Quantitive_Probability_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(prob_plot_QSCAN)
  dev.off()
  
  
  breaks = c()
  for (i in 1:num_true)
  {
    start = Block.ind[locuind[i], 2]; end = Block.ind[locuind[i], 3]
    
    X = start:end; Y = signals[start:end]
    Data = data.frame(X, Y)
    
    cover_plot_dBiRS = cover_plot_dBiRS + geom_point(data = Data, aes(x = X, y = Y), size = 2, col = 'grey') + 
      annotate("rect", xmin = true.ind.frame[i, 1], xmax = true.ind.frame[i, 2], ymin = 0, ymax = Inf, fill= "pink", alpha = 0.4) 
    
    cover_plot_KS = cover_plot_KS + geom_point(data = Data, aes(x = X, y = Y), size = 2, col = 'grey') + 
      annotate("rect", xmin = true.ind.frame[i, 1], xmax = true.ind.frame[i, 2], ymin = 0, ymax = Inf, fill= "pink", alpha = 0.4) 
    
    cover_plot_QSCAN = cover_plot_QSCAN + geom_point(data = Data, aes(x = X, y = Y), size = 2, col = 'grey') + 
      annotate("rect", xmin = true.ind.frame[i, 1], xmax = true.ind.frame[i, 2], ymin = 0, ymax = Inf, fill= "pink", alpha = 0.4) 
    
    breaks = c(breaks, floor((true.ind.frame[i, 1] + true.ind.frame[i, 2])/2))
  }
  
  cover_plot_dBiRS = cover_plot_dBiRS + geom_point(data = Data_dBiRS, aes(x = X_dBiRS, y = Y_dBiRS), size = 2, col = 'red')
  cover_plot_KS = cover_plot_KS + geom_point(data = Data_KS, aes(x = X_KS, y = Y_KS), size = 2, col = 'blue')
  cover_plot_QSCAN = cover_plot_QSCAN + geom_point(data = Data_QSCAN, aes(x = X_QSCAN, y = Y_QSCAN), size = 2, col = 'green')
  
  cover_plot_dBiRS = cover_plot_dBiRS + theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) + 
    scale_y_continuous(breaks = NULL) +
    scale_x_break(c(Block.ind[locuind[1], 3], Block.ind[locuind[2], 2])) +
    scale_x_break(c(Block.ind[locuind[2], 3], Block.ind[locuind[3], 2])) +
    scale_x_break(c(Block.ind[locuind[3], 3], Block.ind[locuind[4], 2])) +
    scale_x_continuous(breaks = breaks) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Underlying Effect")
  
  cover_plot_QSCAN = cover_plot_QSCAN + theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) + 
    scale_y_continuous(breaks = NULL) +
    scale_x_break(c(Block.ind[locuind[1], 3], Block.ind[locuind[2], 2])) +
    scale_x_break(c(Block.ind[locuind[2], 3], Block.ind[locuind[3], 2])) +
    scale_x_break(c(Block.ind[locuind[3], 3], Block.ind[locuind[4], 2])) +
    scale_x_continuous(breaks = breaks) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Underlying Effect")
  
  cover_plot_KS = cover_plot_KS + theme_gray() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) + 
    scale_y_continuous(breaks = NULL) +
    scale_x_break(c(Block.ind[locuind[1], 3], Block.ind[locuind[2], 2])) +
    scale_x_break(c(Block.ind[locuind[2], 3], Block.ind[locuind[3], 2])) +
    scale_x_break(c(Block.ind[locuind[3], 3], Block.ind[locuind[4], 2])) +
    scale_x_continuous(breaks = breaks) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
          axis.title.y = element_text(face = "bold", size = 20), legend.text = element_text(face = "bold", size = 17), 
          axis.text.x = element_text(face = "bold", size = 17),  axis.text.y = element_text(face = "bold", size = 17), legend.position = "none") + 
    xlab("Position") + ylab("Underlying Effect")
  
  png(filename = paste0('Figures/dBiRS_Quantitive_Coverage_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(cover_plot_dBiRS)
  dev.off()
  
  png(filename = paste0('Figures/QSCAN_Quantitive_Coverage_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(cover_plot_QSCAN)
  dev.off()
  
  png(filename = paste0('Figures/KS_Quantitive_Coverage_c = ', cv[setc], '.png'), width = 1200, height = 400)
  print(cover_plot_KS)
  dev.off()
}

