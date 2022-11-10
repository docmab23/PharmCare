arg <- commandArgs(trailingOnly = T)
sdf <- arg[[1]]
pfx <- arg[[2]]
cov <- read.table(sdf, header = F, comment.char = '', colClasses = c('character', 'integer', 'integer'))
pdf(file = paste(pfx, '.pdf', sep = ''))
i <- seq(1, nrow(cov), nrow(cov) / 1000)
plot(cov[i, 2], cov[i, 3], xlab = 'Coordinate', ylab = 'Coverage', ylim = c(0, 100))
lines(lowess(cov[i, 2], cov[i, 3], f = 1/10), col = 'red')
dev.off()