arg <- commandArgs(trailingOnly = T)
dir <- arg[[1]] # program directory
pfx <- arg[[2]] # output prefix
typ <- arg[[3]] # wgs or ts
tgt.nam <- arg[[4]] # target gene
ctl.nam <- arg[[5]] # control gene
ctl.rgn <- arg[[6]] # control region
gdf.nam <- arg[[7]] # GDF file
ctl.typ <- arg[[8]] # control type
ref <- arg[[9]] # reference sample list

# import genes data
gen <- read.table(paste(dir, '/gene_table.txt', sep = ''), header = T, sep = '\t', stringsAsFactors = F)

# get target gene data
tgt.dat <- gen[gen$name == tgt.nam, ]
tgt.chr <- gsub('chr', '', tgt.dat$chr)
tgt.srt <- tgt.dat$hg19_start - tgt.dat$upstream
tgt.end <- tgt.dat$hg19_end + tgt.dat$downstream
tgt.siz <- tgt.end - tgt.srt

# get paralogous gene data
prg.dat <- if (tgt.dat$paralog != '.') gen[gen$name == tgt.dat$paralog, ] else '.'

# get control gene data
ctl.dat <- if (ctl.typ == 'known') gen[gen$name == ctl.nam, ] else '.'
ctl.chr <- gsub('chr', '', gsub(':.*', '', ctl.rgn))
ctl.srt <- as.integer(gsub('-.*', '', gsub('.*:', '', ctl.rgn)))
ctl.end <- as.integer(gsub('.*-', '', gsub('.*:', '', ctl.rgn)))
ctl.siz <- ctl.end - ctl.srt

# read input GDF file
gdf <- read.table(gdf.nam, header = T, check.names = F, stringsAsFactors = F, comment.char = '', sep = '\t', nrows = 5)
gdf <- read.table(gdf.nam, header = T, check.names = F, stringsAsFactors = F, comment.char = '', sep = '\t', colClasses = sapply(gdf, class))

# filter loci with low depth
gdf <- gdf[gdf$Average_Depth_sample >= 5, ]

# divide dataset into two: target gene vs. control gene
names(gdf) <- gsub('Depth_for_', '', names(gdf))
gdf$Total_Depth <- NULL
gdf$Average_Depth_sample <- NULL
chr <- gsub(':.*', '', gsub('chr', '', gdf$Locus))
pos <- as.integer(gsub('.*:', '', gdf$Locus))
tgt <- gdf[chr == tgt.chr & pos >= tgt.srt & pos <= tgt.end, ]
ctl <- gdf[chr == ctl.chr & pos >= ctl.srt & pos <= ctl.end, ]

if (dim(tgt)[1] == 0) {
  stop('sv.R:\n***********************************************************************\nInput GDF file has no data for target gene\n***********************************************************************\n\n')
}

if (dim(ctl)[1] == 0) {
  stop('sv.R:\n***********************************************************************\nInput GDF file has no data for control gene\n***********************************************************************\n\n')
}

tgt$Locus <- as.integer(gsub('.*:', '', tgt$Locus))
ctl$Locus <- as.integer(gsub('.*:', '', ctl$Locus))

# apply intra-sample normalization
cpn <- as.data.frame(cbind(Locus = tgt$Locus, 2 * mapply('/', tgt[, -1, drop = F], colMeans(ctl[, -1, drop = F]))))

# apply inter-sample normalization
if (typ == 'ts') {
  if (ref == '.') { # use population mean
    cpn[, -1] <- 2 * cpn[, -1] / rowMeans(cpn[, -1])
  } else { # use referecne mean
    i <- match(unlist(strsplit(ref, split = ',')), colnames(cpn))
    if (length(i) == 1) { # single sample
      cpn[, -1] <- 2 * cpn[, -1] / cpn[, i]
    } else { # multiple samples
      cpn[, -1] <- 2 * cpn[, -1] / rowMeans(cpn[, i])
    }
  }
}

# divide dataset into two: masked vs. unmasked
if (tgt.dat$masked_starts == '.') {
  msk.dat <- '.'
  msk.pos <- '.'
  ums <- cpn
} else {
  msk.dat <- as.data.frame(cbind(srt = as.integer(unlist(strsplit(tgt.dat$masked_starts, split = ','))), end = as.integer(unlist(strsplit(tgt.dat$masked_ends, split = ',')))))
  msk.pos <- c()
  for (i in 1:nrow(msk.dat)) {
    msk.pos <- c(msk.pos, seq(msk.dat$srt[i], msk.dat$end[i]))
  }
  i <- match(intersect(cpn$Locus, msk.pos), cpn$Locus)
  msk <- cpn[i, ]
  ums <- cpn[-i, ]
}

# apply running medians smoothing
obs <- as.data.frame(cbind(Locus = ums$Locus, sapply(ums[, -1, drop = F], function(x) runmed(x, k = 1001))))

# get SV profiles
svd <- read.table(paste(dir, '/sv_table.txt', sep = ''), header = T, sep = '\t', stringsAsFactors = F)
svd <- svd[svd$gene == tgt.nam, ]
prf <- as.data.frame(tgt.srt:tgt.end)
for (i in 1:length(svd$name)) {
  if (svd$cp[i] == '.') {
    cp <- tgt.siz + 1
  } else {
    cp <- diff(c(tgt.srt - 1, as.numeric(unlist(strsplit(svd$cp[i], split = ','))), tgt.end))
  }
  cn <- as.numeric(unlist(strsplit(svd$cn[i], split = ',')))
  prf <- cbind(prf, rep(cn, cp))
}
colnames(prf) <- append(svd$name, 'Locus', 0)

# make a subsetted copy of SV profiles
i <- match(intersect(prf$Locus, obs$Locus), prf$Locus)
prf.sub <- prf[i, ]

# run svcomb algorithm
seq <- as.data.frame(rbind(t(combn(colnames(prf)[-1], 2)), cbind(colnames(prf)[-1], colnames(prf)[-1])), stringsAsFactors = F)
for (i in 2:ncol(obs)) {
  prob1 <- prob2 <- like <- lambda1 <- lambda2 <- ssr_raw <- ssr <- c(1:nrow(seq))
  for (j in 1:nrow(seq)) { 		
    lambda1[j] <- if (seq[j, 1] != 'no_sv' && seq[j, 2] != 'no_sv') 0.95 else 1
    lambda2[j] <- if (seq[j, 1] == seq[j, 2]) 1 else 0.99
    prob1[j] <- svd$prob[which(svd$name == seq[j, 1])]
    prob2[j] <- svd$prob[which(svd$name == seq[j, 2])]
    like[j] <- prob1[j] * prob2[j]
    ssr_raw[j] <- sum((prf.sub[, seq[j, 1]] + prf.sub[, seq[j, 2]] - obs[, i])^2)
    if (ssr_raw[j] == 0) {
      ssr_raw[j] = 1
    }
    ssr[j] <- round(ssr_raw[j] / like[j] / lambda1[j] / lambda2[j], 1)
    ssr_raw[j] <- round(ssr_raw[j], 1)
  }
  results <- as.data.frame(cbind(seq1 = seq[, 1], seq2 = seq[, 2], prob1, prob2, like, lambda1, lambda2, ssr_raw, ssr))
  write.table(results[order(ssr), ], file = paste(pfx, '.project/ssr/', colnames(obs)[i], '.txt', sep = ''), sep = '\t', quote = F, row.names = F)
}

rm(chr)
rm(pos)
rm(gdf)
rm(cpn)
rm(gen)
rm(msk.pos)
rm(svd)
pfx2 <- pfx
rm(pfx)

# save variables environment
save.image(paste(pfx2, '.project/sv.RData', sep = ''))

warnings()