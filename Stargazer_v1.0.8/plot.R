arg <- commandArgs(trailingOnly = T)
pfx <- arg[[1]]
dtl <- arg[[2]] == 'True'

load(paste(pfx, '.project/sv.RData', sep = ''))

y_max <- 6
y_min <- -1

tgt.avg <- rowMeans(tgt[, -1, drop = F])
ctl.avg <- rowMeans(ctl[, -1, drop = F])

# Get sample data
result <- read.table(paste(pfx, '.txt', sep = ''), header = T, as.is = T, sep = '\t', comment.char = '', colClasses = c('character'))

# set x-axis increments
stp <- if (dtl) 1000 else 200
tgt.cov <- nrow(tgt) / tgt.siz
ctl.cov <- nrow(ctl) / ctl.siz
tgt.inc <- seq(1, tgt.siz, round(tgt.siz / stp * tgt.cov))
ctl.inc <- seq(1, ctl.siz, round(ctl.siz / stp * ctl.cov))

# define plotting functions
plot_base <- function(r, ytitle) {
  plot(x = 0, type = 'n', xlim = c(tgt.srt, tgt.end), ylim = c(y_min * r, y_max * r), xlab = paste('Chromosome ', tgt.chr, sep =''),	ylab = '', yaxt = 'n')
  title(ylab = ytitle)
  add_gene_model(tgt.dat, r)
  add_gene_name(tgt.dat$hg19_start, tgt.dat$hg19_end, tgt.nam, r)
  if (tgt.dat$paralog != '.') {
    add_gene_model(prg.dat, r)
    add_gene_name(prg.dat$hg19_start, prg.dat$hg19_end, prg.dat$name, r)
  }
}

add_gene_model <- function(g, r) {
	lines(x = c(g$hg19_start, g$hg19_end), y = c(y_min * 0.5 * r, y_min * 0.5 * r))
	rect(xleft = as.numeric(unlist(strsplit(g$exon_starts, split = ','))), ybottom = y_min * 0.5 * r + 0.15 * r, xright = as.numeric(unlist(strsplit(g$exon_ends, split = ','))), ytop = y_min * 0.5 * r - 0.15 * r, col = 'black', border = 'black')
}

add_gene_name <- function(start, end, name, scale) {
	text(x = (start + end) / 2, y = -0.9 * scale, labels = bquote(italic(.(toupper(name)))), cex = 0.83)
}

add_masked_regions <- function(r) {
  if (is.character(msk.dat)) return() else rect(xleft = msk.dat[, 1], ybottom = 0, xright = msk.dat[, 2], ytop = y_max * r, col = rgb(0, 0, 0, 0.25), border = NA)
}

add_svcomb <- function(sc1, sc2) {
	points(prf[, 1], prf[, sc1] + prf[, sc2], type = 'l', col = 'red', lwd = 6)
	points(prf[, 1], prf[, sc1], type = 'l', col = 'navy', lwd = 2)
	points(prf[, 1], prf[, sc2], type = 'l', col = 'cyan', lwd = 2, lty = 2)
}

# make plots for every sample
for (i in 1:nrow(result)) {
	id <- result[i, 'name']
	hap1_sv <- result[i, 'hap1_sv']
	hap2_sv <- result[i, 'hap2_sv']

	if (hap1_sv == '.') {
		dip_sv <- unlist(strsplit(result[i, 'dip_sv'], split = ','))
		hap1_sv <- dip_sv[1]
		hap2_sv <- dip_sv[2]
	}	
	
	# open graphics device
	pdf(file = paste(pfx, '.project/plot/', id, '.pdf', sep = ''))
	par(mfrow = c(2, 2))

	# determine y-axis scale for read depth profiles
	sc1 <- if (mean(tgt[, i + 1]) < 60) 15 else 30

	# plot read depth profile for target gene
	plot_base(sc1, 'Read depth')
	axis(2, at = seq(0, y_max * sc1, 2 * sc1))
	points(tgt[tgt.inc, 1], tgt[tgt.inc, i + 1], pch = 21, col = 'black', bg = 'gray')
	points(tgt[tgt.inc, 1], tgt.avg[tgt.inc], pch = 21, col = 'black', bg = 'green')

	# plot read depth profile for control gene
	plot(x = 0, type = 'n', xlim = c(ctl.srt, ctl.end), ylim = c(y_min * sc1, y_max * sc1), xlab = paste('Chromosome ', ctl.chr, sep = ''),	ylab = '', yaxt = 'n')
	title(ylab = 'Read depth')
	axis(2, at = seq(0, y_max * sc1, 2 * sc1))
	points(ctl[ctl.inc, 1], ctl[ctl.inc, i + 1], pch = 21, col = 'black', bg = 'gray')
	points(ctl[ctl.inc, 1], ctl.avg[ctl.inc], pch = 21, col = 'black', bg = 'green')
	if (ctl.typ == 'known') {
		add_gene_model(ctl.dat, sc1)
		add_gene_name(ctl.srt, ctl.end, ctl.nam, sc1)	
	} else {
		add_gene_name(ctl.srt, ctl.end, ctl.typ, sc1)
	}

	# determine y-axis scale for copy number profile and allele fraction profile
	sc2 <- if (max(prf[, hap1_sv] + prf[, hap2_sv]) < 5) 1 else 2

	# plot copy number profile 
	plot_base(sc2, 'Copy number')
	axis(2, at = seq(0, y_max * sc2, 2 * sc2))
	if (nrow(ums) == 0) {
		return()
	}
	add_svcomb(hap1_sv, hap2_sv)
	if (!is.character(msk.dat)) {
	  points(msk[tgt.inc, 1], msk[tgt.inc, i + 1], pch = 21, col = 'black', bg = 'gray')
	}
	points(ums[tgt.inc, 1], ums[tgt.inc, i + 1], pch = 21, col = 'black', bg = 'gray')
	if (dtl) {
		points(obs[tgt.inc, 1], obs[tgt.inc, i + 1], pch = 21, col = 'black', bg = 'deepskyblue')
		add_masked_regions(sc2)
	}

	# plot allele fraction profile
	af <- read.table(paste(pfx, '.project/af/', id, '.txt', sep = ''), header = T)
	plot_base(sc2, 'Allele fraction')
	axis(2, at = c(0, y_max * sc2 / 2, y_max * sc2), labels = c(0, 0.5, 1))
	axis(4, at = seq(0, y_max * sc2, 2 * sc2))
	add_svcomb(hap1_sv, hap2_sv)
	if (dim(af)[1] == 0) {
		dev.off()
		next
	}
	points(af[, 1], af[, 'hap1_af'] * y_max * sc2, pch = 21, col = 'black', bg = 'navy')
	points(af[, 1], af[, 'hap2_af'] * y_max * sc2, pch = 21, col = 'black', bg = 'cyan')
	if (dtl) {
		add_masked_regions(sc2)
	}

	# shut down grpahics device
	dev.off()
}

warnings()

