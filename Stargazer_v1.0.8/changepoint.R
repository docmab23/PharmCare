library(changepoint)

arg <- commandArgs(trailingOnly = T)
pfx <- arg[[1]]
svd <- arg[[2]]

load(svd)

cpt <- cpt.meanvar(
  data = as.matrix(t(obs[, -1, drop = F])),
  method = 'BinSeg',
  test.stat = 'Exponential',
  minseglen = 500
)

cat(paste(paste('sample', 'changepoints', 'means', sep = '\t'), '\n'), file = paste(pfx, '.txt', sep = ''))

for (i in 1:length(cpt)) {
  sid <- names(obs)[i + 1] # sample ID
  pts <- paste(paste(obs[cpt[[i]]@cpts, 1], collapse = ','), ',', sep = '') # detected changepoints
  avg <- paste(paste(round(1 / cpt[[i]]@param.est$rate, 2), collapse = ','), ',', sep = '') # detected means
  line <- paste(sid, pts, avg, sep = '\t')
  cat(paste(line, '\n'), file = paste(pfx, '.txt', sep = ''), append = T)
}

warnings()

