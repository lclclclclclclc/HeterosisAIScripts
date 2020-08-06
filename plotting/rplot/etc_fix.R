setwd("~/Documents/science/Grad/demog20/proj/HeterosisAIScripts")
title <- c("n","insert_ai", "growth","meanpI","start","end","freq_before","freq_after","pI","Dstat","fD","Het","divratioavg","Q_1_100_q95","Q_1_100_q90","Q_1_100_max","U_1_0_100","U_1_20_100","U_1_50_100","U_1_80_100")
sexX_null <- read.table("output/stats/sexXnull_nscale100/sgcz-dominance0-model0-sexX-hs0-ai0attempt645forItAll_human_windows.txt")
sexX_ai <- read.table("output/stats/sexX_nscale100/sgcz-dominance0-model0-sexX-hs0-ai0.01_human_windows.txt")

colnames(sexX_null) <- title
colnames(sexX_ai) <- title


sexX_null$pos <- round((sexX_null$end+sexX_null$start)/2)
sexX_ai$pos <- round((sexX_ai$end+sexX_ai$start)/2)


plot(sexX_ai$pos, sexX_ai$fD)
plot(sexX_null$pos, sexX_null$fD)

plot(sexX_ai$pos, sexX_ai$U_1_80_100)
plot(sexX_null$pos, sexX_null$U_1_80_100)

