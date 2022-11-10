import copy
import gzip
import os
import shutil
import statistics
import errno

PHENOTYPES = ['unfavorable_response', 'poor_function', 'poor_metabolizer', 'decreased_function', 'slow_metabolizer', 'intermediate_metabolizer', 'normal_function', 'normal_metabolizer', 'rapid_metabolizer', 'increased_function', 'ultrarapid_metabolizer', 'favorable_response', 'unknown', '.']

class VCF:
	def __init__(self):
		self.meta = []
		self.header = []
		self.data = []

class GDF:
	def __init__(self):
		self.header = []
		self.data = []

class SNP:
	def __init__(self):
		self.pos = '' # reference genome position
		self.wt = '' # wild type (*1) allele
		self.var = '' # variant allele
		self.rs = '' # rs ID
		self.het = False # heterozygous
		self.ad = 0 # allelic depth
		self.td = 0 # total depth
		self.n = '' # SNP table number
		self.hg = '' # reference genome allele
		self.so = '' # sequence ontology
		self.effect = '' # coding effect
		self.impact = '' # variant impact
		self.rev = False # reverting variant

	@property
	def key(self):
		return (self.pos, self.wt, self.var)

	@property
	def af(self):
		return 0 if self.td == 0 else self.ad / self.td
	
	def __eq__(self, other):
		return self.key == other.key

	def __hash__(self):
		return hash(self.key)
		
	def summary(self):	
		return '<{}:{}>{}:{}/{}:{:.2f}:{}:{}:{}>'.format(self.pos, self.wt, self.var, self.ad, self.td, self.af, self.so, self.impact, self.effect)

class Star:
	def __init__(self):
		self.name = ''
		self.score = -1.0
		self.core = []
		self.tag = []
		self.sv = ''

	@property
	def ranked_as(self):
		if self.score == "unknown": # when sorted, unknown function alleles should be broken ties with normal function alleles using other attributes
			return 1.0
		elif self.score > 1: # when sorted, increased function alleles should come before normal function alleles
			return 0.99
		else:
			return self.score

	@property
	def rank(self):
		return (self.ranked_as, -1 * int(bool(self.sv)), -1 * len(self.core))

	def __str__(self):
		return self.name

	def __repr__(self):
		return str(self)

	def __eq__(self, other):
		return self.name == other.name

	def __hash__(self):
		return hash(self.name)

class Haplotype:
	def __init__(self):
		self.cand = []
		self.obs = []

	@property
	def sv(self):
		sv = 'no_sv'
		sv_list = []
		for star_allele in self.cand:
			if star_allele.sv and star_allele.sv not in sv_list:
				sv_list.append(star_allele.sv)
		if len(sv_list) > 1:
			raise ValueError('haplotype contains multiple structural variant calls')
		if len(sv_list) == 1:
			sv = sv_list[0]
		return sv

	@property
	def af(self):
		return [0] if not self.obs else [x.af for x in self.obs]

	@property
	def af_mean_main(self):
		filtered = [x.af for x in self.obs if x.td > 10 and x.het and x in [y for y in self.cand[0].core]]
		return -1 if not filtered else statistics.mean(filtered)

	def af_mean_gene(self, start, end):
		filtered = [x.af for x in self.obs if x.td > 10 and x.het and start <= int(x.pos) <= end]
		return -1 if not filtered else statistics.mean(filtered)
		
	def fit_data(self, total_cn, start, end):
		"""Return the fit MAF and CN."""
		maf_choices = []
		for i in range(1, total_cn):
			maf_choices.append(i / total_cn)
		fit_maf = min(maf_choices, key = lambda x: abs(x - self.af_mean_gene(start, end)))
		fit_cn = maf_choices.index(fit_maf) + 1
		return fit_maf, fit_cn

	def remove_star(self, sX):
		"""Remove the given star allele from the candidates list."""
		for i, star in enumerate(self.cand):
			if star.name == sX.name:
				del self.cand[i]
				break

	def add_dup(self, cn):
		"""Duplicate the main star allele by the given CN."""
		if cn == 1:
			return
		if cn > 10:
			cn = 10
		sX = self.cand[0]
		if isinstance(sX.score, str):
			score = 'unknown'
		else:
			score = sX.score * cn
		name = sX.name + 'x' + str(cn)
		sY = Star()
		
		sY.name = name; sY.score = score; sY.core = copy.deepcopy(sX.core); sY.sv = 'cnv{}'.format(cn)
		
		self.cand.insert(0, sY)
		self.remove_star(sX)

class Sample:
	def __init__(self):
		self.name = '' # sample ID
		self.gt = False # true if genotyped
		self.sv = ['', ''] # SV calls
		self.pt = '' # predicted phenotype
		self.ssr = '' # sum of squared residuals
		self.dip_cand = [] # candidate stars
		self.hap = [Haplotype(), Haplotype()]
		self.bad = False # true if QC failed

	@property
	def dip_score(self):
		hap_scores = [self.hap[0].cand[0].score, self.hap[1].cand[0].score]	
		if 'unknown' in hap_scores:
			return 'unknown'
		else:
			return sum(hap_scores)

def copy_vcf(original_vcf, items):
	copied_vcf = VCF()
	if 'meta' in items:
		copied_vcf.meta = copy.deepcopy(original_vcf.meta)
	if 'header' in items:
		copied_vcf.header = copy.deepcopy(original_vcf.header)
	if 'data' in items:
		copied_vcf.data = copy.deepcopy(original_vcf.data)
	return copied_vcf

def parse_region(region):
	return {'chr': region.split(':')[0].replace('chr', ''), 'start': int(region.split(':')[1].split('-')[0]), 'end': int(region.split(':')[1].split('-')[1])}

def read_vcf_simple(file):
	f = gzip.open(file, 'rt') if '.gz' in file else open(file)
	vcf = VCF()
	for line in f:
		if '##' in line:
			vcf.meta.append(line)
			continue
		fields = line.strip().split('\t')
		if fields[0] == '#CHROM':
			vcf.header = fields
			continue
		chr = fields[0].replace('chr', '')
		vcf.data.append([chr] + fields[1:])
	f.close()
	return vcf

def read_vcf_region(file, region):
	vcf = VCF()
	region_dict = parse_region(region)
	f = gzip.open(file, 'rt') if '.gz' in file else open(file)
	for line in f:
		if '##' in line:
			vcf.meta.append(line)
			continue
		fields = line.strip().split('\t')
		if fields[0] == '#CHROM':
			vcf.header = fields
			continue
		chr, pos = fields[0].replace('chr', ''), int(fields[1])
		if chr != region_dict['chr'] or pos < region_dict['start']:
			continue
		if pos > region_dict['end']:
			break
		vcf.data.append([chr] + fields[1:])
	f.close()
	return vcf

def write_vcf(vcf, file):
	with open(file, 'w') as f:
		for line in vcf.meta:
			f.write(line)
		f.write('\t'.join(vcf.header) + '\n')
		for fields in vcf.data:
			f.write('\t'.join(fields) + '\n')

def write_gdf(gdf, file):
	with open(file, 'w') as f:
		f.write('\t'.join(gdf.header) + '\n')
		for fields in gdf.data:
			f.write('\t'.join([str(x) for x in fields]) + '\n')

def get_gene_dict():
	gene_dict = {}
	with open(os.path.dirname(os.path.realpath(__file__)) + '/gene_table.txt') as f:
		header = next(f).strip().split('\t')
		for line in f:
			fields = line.strip().split('\t')
			name = fields[header.index('name')]
			chr = fields[header.index('chr')].replace('chr', '')
			hg19_start = int(fields[header.index('hg19_start')])
			hg19_end = int(fields[header.index('hg19_end')])
			upstream = int(fields[header.index('upstream')])
			downstream = int(fields[header.index('downstream')])
			region = '{}:{}-{}'.format(chr, hg19_start - upstream, hg19_end + downstream)
			gene_dict[name] = dict(zip(header + ['region'], fields + [region]))
	return gene_dict

def get_snp_list(target_gene):
	snp_list = []
	with open(os.path.dirname(os.path.realpath(__file__)) + '/snp_table.txt') as f:
		next(f)
		for line in f:
			fields = line.strip().split('\t')
			gene = fields[0]
			if gene != target_gene:
				continue
			snp = SNP()
			snp.n, snp.effect, snp.pos, snp.id, snp.hg, snp.var, snp.wt, snp.so, snp.impact, snp.rev = fields[1], fields[4], fields[5], fields[6], fields[7], fields[8], fields[9], fields[10], fields[11], fields[12] == 'yes'
			snp_list.append(snp)
	return snp_list

def get_star_dict(target_gene, snp_list):
	star_dict = {}
	with open(os.path.dirname(os.path.realpath(__file__)) + '/star_table.txt') as f:
		next(f)
		for line in f:
			fields = line.strip().split('\t')
			gene = fields[0]
			if gene != target_gene:
				continue
			star = Star()
			star.name = fields[2]
			star.score = float(fields[8]) if fields[8] != 'unknown' else fields[8]
			star.core = [] if fields[4] == 'ref' or fields[4] == '.' else copy.deepcopy([x for x in snp_list if '{}:{}>{}'.format(x.pos, x.wt, x.var) in fields[4].split(',')])
			star.tag = [] if fields[5] == '.' else copy.deepcopy([x for x in snp_list if '{}:{}>{}'.format(x.pos, x.wt, x.var) in fields[5].split(',')])
			star.sv = '' if fields[6] == '.' else fields[6]
			star_dict[star.name] = star
	return star_dict

def vcf2samples(vcf):
	samples = {}
	for name in vcf.header[9:]:
		sample = Sample()
		sample.name = name
		i = vcf.header.index(name)
		for fields in vcf.data:
			pos, rs, ref, alt, inf, fmt = fields[1], fields[2], fields[3], fields[4].split(','), fields[7].split(';'), fields[8]
		
			if not any(['PS=D' in x for x in inf]):
				continue

			gt = [int(x) for x in fields[i].split(":")[0].split("|")]
			al = [ref] + alt
			
			vi_list = ['no_change'] + [x for x in inf if 'VI=' in x][0].replace('VI=', '').split(',')
			so_list = ['no_change'] + [x for x in inf if 'SO=' in x][0].replace('SO=', '').split(',')
			fe_list = ['no_change'] + [x for x in inf if 'FE=' in x][0].replace('FE=', '').split(',')
			
			for j in [0, 1]:
				
				snp = SNP()
				snp.pos, snp.wt, snp.var, snp.rs, snp.het, snp.so, snp.impact, snp.effect = pos, ref, al[gt[j]], rs, gt[0] != gt[1], so_list[gt[j]], vi_list[gt[j]], fe_list[gt[j]]
				
				
				if 'AD' in fmt:
					ad_list = [int(x) for x in fields[i].split(':')[1].split(',')]
					snp.ad = ad_list[gt[j]]; snp.td = sum(ad_list)
					
				sample.hap[j].obs.append(snp)
		samples[name] = sample		
	return samples

def delete_file(file):
	try:
		os.remove(file)
	except OSError:
		pass

def delete_dir(dir):
	try:
		shutil.rmtree(dir)
	except OSError:
		pass

def create_dir(dir):
	delete_dir(dir)
	os.mkdir(dir)

def sort_star_names(names):
	def f(x):
		cn = 1
		if '*' not in x or x == '*DEL':
			n = 999
		else:
			n = int(''.join([y for y in x.split('+')[0].split('x')[0] if y.isdigit()]))
			if 'x' in x.split('+')[0]:
				cn = int(x.split('+')[0].split('x')[1])
		return (n, cn, len(x))

	return sorted(names, key = f)

def check_file(file):
	if os.path.isfile(file) or os.path.isdir(file):
		return os.path.realpath(file)
	raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def sample_regions(region):
	size = 1000
	region_dict = parse_region(region)
	center = int((region_dict['start'] + region_dict['end']) / 2)
	start = '{}:{}-{}'.format(region_dict['chr'], region_dict['start'], region_dict['start'] + size)
	middle = '{}:{}-{}'.format(region_dict['chr'], round(center - size / 2), round(center + size / 2))
	end = '{}:{}-{}'.format(region_dict['chr'], region_dict['end'] - size, region_dict['end'])
	return [start, middle, end]
	
def sort_regions(regions):
	def f(x):
		region_dict = parse_region(x)
		chr = 23 if region_dict['chr'] == 'X' else 24 if region_dict['chr'] == 'Y' else int(region_dict['chr'])
		return (chr, region_dict['start'], region_dict['end'])
	return sorted(regions, key = f)
	
def get_target_genes():
	return [x for x in get_gene_dict() if get_gene_dict()[x]['type'] == 'target']

def get_control_genes():
	return [x for x in get_gene_dict() if get_gene_dict()[x]['control'] == 'yes']
	
def get_cnsr(gene):
	return get_gene_dict()[gene]['cnsr'].replace('chr', '')
	
def get_region(gene):
	return get_gene_dict()[gene]['region']

def get_paralog(gene):
	return get_gene_dict()[gene]['paralog']
	
def get_function(gene):
	return get_gene_dict()[gene]['function'].capitalize()

def get_pv_member(gene):
	return get_gene_dict()[gene]['pv_member'].capitalize()

def get_dpsv_member(gene):
	return get_gene_dict()[gene]['dpsv_member'].capitalize()
	
def get_chr(gene):
	return get_gene_dict()[gene]['chr'].replace('chr', '')
	
def get_masked_starts(gene):
	return [int(x) for x in get_gene_dict()[gene]['masked_starts'].split(',')[:-1]]

def get_masked_ends(gene):
	return [int(x) for x in get_gene_dict()[gene]['masked_ends'].split(',')[:-1]]
	
def get_hg19_start(gene):
	return int(get_gene_dict()[gene]['hg19_start'])
	
def get_hg19_end(gene):
	return int(get_gene_dict()[gene]['hg19_end'])