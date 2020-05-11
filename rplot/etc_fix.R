setwd("~/Documents/science/Grad/demog20/proj/HeterosisAIScripts")
title <- c("n","insert_ai", "growth","meanpI","start","end","freq_before","freq_after","pI","Dstat","fD","Het","divratioavg","Q_1_100_q95","Q_1_100_q90","Q_1_100_max","U_1_0_100","U_1_20_100","U_1_50_100","U_1_80_100")
neu <- read.table("output/stats/chr19region-dominance0.01-model0-growth4-hs0-ai0.0_human_windows.txt")
neg <- read.table("output/stats/chr19region-dominance2-model0-growth4-hs0-ai0.0_human_windows.txt")
neg_sel <- read.table("output/stats/chr19region-dominance2-model0-growth4-hs0-ai0.1_human_windows.txt")

colnames(neu) <- title
colnames(neg) <- title
colnames(neg_sel) <- title
neg$window <- c(1:99)
neg$pos <- round((neg$end+neg$start)/2)
neu$window <- c(1:99)
neu$pos <- round((neu$end+neu$start)/2)
neg_sel$window <- c(1:99)
neg_sel$pos <- round((neg_sel$end+neg_sel$start)/2)
