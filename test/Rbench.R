library(vars)
library(frequencyConnectedness)

i <- as.integer(commandArgs(trailingOnly = TRUE)[1])
H <- as.integer(commandArgs(trailingOnly = TRUE)[2])
w <- as.integer(commandArgs(trailingOnly = TRUE)[3])
l <- as.integer(commandArgs(trailingOnly = TRUE)[4])


data <- read.csv("data.csv", header = F)
data <- log(data)

bounds <- c(pi + 0.001, pi/5, pi/10, 0)

est <- VAR(data[1:w + i,], p = l, type = "const")

write.table(do.call(rbind,spilloverBK12(est, n.ahead = H, no.corr = T, partition = bounds, table = T, absolute = T)), file = "table_T.txt", sep = ",", row.names = F, col.names = F)
write.table(do.call(rbind,spilloverBK12(est, n.ahead = H, no.corr = F, partition = bounds, table = T, absolute = T)), file = "table_F.txt", sep = ",", row.names = F, col.names = F)
write.table(spilloverBK12(est, n.ahead = H, no.corr = F, partition = bounds, table = F, absolute = T), file = "overall_FT.txt", sep = ",", row.names = F, col.names = F)
write.table(spilloverBK12(est, n.ahead = H, no.corr = T, partition = bounds, table = F, absolute = T), file = "overall_TT.txt", sep = ",", row.names = F, col.names = F)
write.table(spilloverBK12(est, n.ahead = H, no.corr = F, partition = bounds, table = F, absolute = F), file = "overall_FF.txt", sep = ",", row.names = F, col.names = F)
write.table(spilloverBK12(est, n.ahead = H, no.corr = T, partition = bounds, table = F, absolute = F), file = "overall_TF.txt", sep = ",", row.names = F, col.names = F)