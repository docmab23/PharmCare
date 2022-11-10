import copy
import helper
import subprocess
import statistics
import operator
import os

DESCRIPTION = f'''tool description:
  call star alleles in target gene from genomic data

getting help:
  stargazer.py genotype -h

main usages:
  genotype with WGS data (or TS data)
    stargazer.py genotype -o OUTPUT_PREFIX -d wgs (or -d ts) -t TARGET_GENE --vcf VCF -c CONTROL_GENE --gdf GDF

  genotype with WGS-generated (or TS-generated) VCF only
    stargazer.py genotype -o OUTPUT_PREFIX -d wgs (or -d ts) -t TARGET_GENE --vcf VCF

  genotype with SNP array-generated VCF only
    stargazer.py genotype -o OUTPUT_PREFIX -d chip -t TARGET_GENE --vcf VCF

other usages:
  genotype with WGS data (or TS data), using the CNSR as the control locus
    stargazer.py genotype -o OUTPUT_PREFIX -d wgs (or -d ts) -t TARGET_GENE --vcf VCF --gdf GDF --control_type cnsr

  genotype with WGS data (or TS data), using a custom region as the control locus
    stargazer.py genotype -o OUTPUT_PREFIX -d wgs (or -d ts) -t TARGET_GENE --vcf VCF --gdf GDF --control_type custom --region REGION

  genotype with TS data, using one or more known SV-free samples as the control set for inter-sample normalization of read depth
    stargazer.py genotype -o OUTPUT_PREFIX -d ts -t TARGET_GENE --vcf VCF -c CONTROL_GENE --gdf GDF --sample_list [SAMPLE [SAMPLE ...]]

  genotype with SNP array-generated VCF only, using the Beagle program's imputation algorithm
    stargazer.py genotype -o OUTPUT_PREFIX -d chip -t TARGET_GENE --vcf VCF --impute
'''

def cyp2a6(sample):
	call_sv1(sample, "gc_e1e4", "*34")
	call_sv1(sample, "gc_e1e2", "*12")
	call_tandem(sample, "dup2", "*2", "*S1")
	call_tandem(sample, "dup2", "*1", "*S1")
	call_tandem(sample, "gc_e9", "*1", "*S2")
	call_tandem(sample, "dup7", "*1", "*S3")
	call_tandem(sample, "dup7b", "*1", "*S6")
	cyp2a6_svcomb(sample)

def cyp2b6(sample):
	call_sv1(sample, "gc_i4e9", "*29")

def cyp2d6(sample):
	call_tandem(sample, 'del1', '*S2', '*1', ordered = True)
	call_tandem(sample, 'gc_i1e9_5', '*68x5', '*4')
	call_tandem(sample, 'gc_i1e9', '*68', '*4')
	call_tandem(sample, 'gc_i1e9', '*S1', '*1')
	call_tandem(sample, 'gc_e9', '*4N', '*4')
	call_tandem(sample, 'gc_e9_3', '*36x3', '*10')
	call_tandem(sample, 'gc_e9_7', '*36x7', '*10')
	call_tandem(sample, 'gc_e9', '*36', '*10')
	call_tandem(sample, 'gc_e9', '*83', '*2')
	call_tandem(sample, 'gc_7to6_i4', '*13A', '*2')
	call_sv1(sample, 'gc_e1e7', '*13C')
	call_sv1(sample, 'gc_7to6_i1', '*13B')
	cyp2d6_svcomb(sample)

def cyp2e1(sample):
	call_sv1(sample, "dup_e7e9", "*S1")

def gstm1(sample):
	gstm1_svcomb(sample)

def gstt1(sample):
	gstt1_svcomb(sample)

def slc22a2(sample):
	call_sv1(sample, "del_i9", "*S1")
	call_sv1(sample, "del_e11", "*S2")

def slco1b1(sample):
	call_sv1(sample, "dup1", "*S3")

def ugt1a4(sample):
	call_sv1(sample, 'del_i1', '*S1')
	call_sv1(sample, 'del2', '*S2')

def ugt2b15(sample):
	call_sv1(sample, "del_i3e6", "*S1")

def ugt2b17(sample):
	ugt2b17_svcomb(sample)

def new_tandem(sv, star_names):
	stars = [args_.star_dict[x] for x in star_names]
	score = sum([x.score for x in stars]) if all([isinstance(x.score, float) for x in stars]) else 'unknown'
	core = list(set([x for y in [z.core for z in stars] for x in y]))
	tandem = helper.Star()
	tandem.name, tandem.score, tandem.core, tandem.sv = '+'.join(star_names), score, copy.deepcopy(core), sv
	return tandem
	
def new_dup(sX, cnv):
	times = int(cnv.replace('cnv', ''))		
	score = 'unknown' if isinstance(sX.score, str) else sX.score * times
	dup = helper.Star()	
	dup.name, dup.score, dup.core, dup.sv = sX.name + 'x' + str(times), score, copy.deepcopy(sX.core), cnv
	return dup

def remove_select(hap, stars):
	'''
	args_:
		hap (list of Stars)
		stars (list of Stars)
	'''
	for i in reversed(range(len(hap))):
		if hap[i].name in [x.name for x in stars]:
			del hap[i]

def remove_sv(hap, l = []):
	for i in reversed(range(len(hap))):
		if hap[i].sv:
			if hap[i].name in l:
				continue
			else:
				del hap[i]
			
def which_has(sample, stars):
	'''
	args_:
		sample (Sample)
		stars (list of str)
	Returns:
		i (int)
	'''
	h1 = set(stars).issubset([x.name for x in sample.hap[0].cand])		
	h2 = set(stars).issubset([x.name for x in sample.hap[1].cand])
	if h1 and h2:
		i = 3
	elif h1 and not h2:
		i = 1
	elif not h1 and h2:
		i = 2
	else:
		i = 0
	return i

def which_severe(sample):
	'''
	args_:
		sample (Sample)
	Returns:
		i (int)
	'''
	h1 = copy.deepcopy(sample.hap[0].cand)
	h2 = copy.deepcopy(sample.hap[1].cand)
	remove_sv(h1)
	remove_sv(h2)
	if h1[0].name == h2[0].name:
		i = 3
	else:
		h3 = sorted([h1[0], h2[0]], key = lambda x: x.rank)
		if h3[0].name == h1[0].name:
			i = 1
		else:
			i = 2
	return i

def call_sv1(sample, sv, x):
	'''
	This function calls a sample's final genotype if the sample has only one SV.
	If a SV-carrying allele is in LD with other alleles, the funtion takes those other alleles as input. 
	x = the name of the star allele with the SV
	'''
	if sample.gt or sample.sv != ["no_sv", sv]:
		return
	if args_.star_dict[x].core:		
		h1 = set(args_.star_dict[x].core).issubset(sample.hap[0].obs)
		h2 = set(args_.star_dict[x].core).issubset(sample.hap[1].obs)
	else:
		h1 = True
		h2 = True
	if not h1 and not h2:
		return
	elif h1 and not h2:
		i, j = 1, 0
	elif not h1 and h2:
		i, j = 0, 1
	else:		
		l1 = copy.deepcopy(sample.hap[0].cand)
		l2 = copy.deepcopy(sample.hap[1].cand)
		remove_sv(l1)
		remove_sv(l2)
		if l1[0].name == l2[0].name:
			i, j = 0, 1
		else:
			l3 = sorted([l1[0], l2[0]], key = lambda x: x.rank)
			if l3[0].name == l1[0].name:
				i, j = 0, 1
			else:
				i, j = 1, 0
	remove_sv(sample.hap[i].cand)
	remove_sv(sample.hap[j].cand, [x])
	sample.gt = True

def call_tandem(sample, sv, x, y, ordered = False):
	"""
	Calls a tandem duplication allele containing two gene copies (e.g., CYP2D6*36+*10)
	x = the name of the 1st star allele in the tandem (e.g., '*36')
	y = the name of the 2nd star allele in the tandem (e.g., '*10')
	"""
	if sample.gt or sample.sv != ["no_sv", sv]:
		return		
	h1 = set([x, y]).issubset([_.name for _ in sample.hap[0].cand])
	h2 = set([x, y]).issubset([_.name for _ in sample.hap[1].cand])
	if not h1 and not h2:
		return
	elif h1 and not h2:
		i, j = 0, 1
	elif not h1 and h2:
		i, j = 1, 0
	else:
		l1 = copy.deepcopy(sample.hap[0].cand)
		l2 = copy.deepcopy(sample.hap[1].cand)
		remove_sv(l1)			
		remove_sv(l2)				
		if l1[0] == l2[0]:
			i, j = 1, 0
		else:
			l3 = sorted([l1[0], l2[0]], key = lambda x: x.rank)
			if l3[0] == l1[0]:
				i, j = 1, 0	
			else:
				i, j = 0, 1
	sX = args_.star_dict[x]
	sY = args_.star_dict[y]
	
	# find SNPs shared by both star alleles
	overlap = []
	for snp in sX.core:
		if snp in sY.core:
			overlap.append(snp)
	
	# return if allele fraction in any of the shared SNPs is less than 0.4
	for snp in overlap:
		if [x for x in sample.hap[i].obs if x == snp][0].af < 0.4:
			return	
	
	tandem = new_tandem(sv, [sX.name, sY.name])			
	sample.hap[i].cand.insert(0, tandem)
	remove_select(sample.hap[i].cand, [sX, sY])		
	remove_sv(sample.hap[i].cand, [tandem.name])
	remove_sv(sample.hap[j].cand)
	if ordered:
		sample.hap[i].cand.sort(key = lambda x: x.rank)
	sample.gt = True

def call_cnv3(sample):
	if sample.gt or sample.sv != ["no_sv", 'cnv2']:
		return
	sX = sample.hap[0].cand[0]
	sY = sample.hap[1].cand[0]

	# Simplest case (e.g., *2/*2x2)
	if sX.name == sY.name:
		sample.hap[0].add_dup(2)
		sample.gt = True
		return

	sX_gene = sample.hap[0].af_mean_gene(helper.get_hg19_start(args_.target_gene), helper.get_hg19_end(args_.target_gene))
	sY_gene = sample.hap[1].af_mean_gene(helper.get_hg19_start(args_.target_gene), helper.get_hg19_end(args_.target_gene))
	if sX_gene == -1:
		sX_gene = 1 - sY_gene
	if sY_gene == -1:
		sY_gene = 1 - sX_gene
	diff_gene = sY_gene - sX_gene

	sX_main = sample.hap[0].af_mean_main
	sY_main = sample.hap[1].af_mean_main		
	if sX_main == -1:
		 sX_main = 1 - sY_main
	if sY_main == -1:
		sY_main = 1 - sX_main
	diff_main = sY_main - sX_main	
	
	
	fit_maf1, _ = sample.hap[0].fit_data(3, helper.get_hg19_start(args_.target_gene), helper.get_hg19_end(args_.target_gene))
	fit_maf2, _ = sample.hap[1].fit_data(3, helper.get_hg19_start(args_.target_gene), helper.get_hg19_end(args_.target_gene))
	
		
	means = [round(fit_maf1, 2), round(fit_maf2, 2)]


	f = lambda a, b: ((a == b) & (a == 0)) | (a * b > 0)
		
	if f(diff_gene, diff_main):
		if means == [0.33, 0.67]:
			sample.hap[1].add_dup(2)
			sample.gt = True
		elif means == [0.67, 0.33]:
			sample.hap[0].add_dup(2)
			sample.gt = True
	else:
		
		if abs(diff_main) > abs(diff_gene):
			
			if sY_main > sX_main:
				sample.hap[1].add_dup(2)
				sample.gt = True
				
			else:
				sample.hap[0].add_dup(2)
				sample.gt = True
		
		else:
			if means == [0.33, 0.67]:
				sample.hap[1].add_dup(2)
				sample.gt = True
			elif means == [0.67, 0.33]:
				sample.hap[0].add_dup(2)
				sample.gt = True		

def call_cnv_plus(sample):
	'''This function calls a final genotype with CN > 3 gene copies.'''
	if sample.gt:
		return
	if 'cnv' not in sample.sv[1]:
		return
	if sample.sv[0] == 'no_sv' and (sample.sv[1] == 'cnv0' or sample.sv[1] == 'cnv2'):
		return
	if sample.sv[0] == 'no_sv' and 'cnv' in sample.sv[1]:
		total_cn = int(sample.sv[1].replace('cnv', '')) + 1
	elif 'cnv' in sample.sv[0] and 'cnv' in sample.sv[1]:
		total_cn = int(sample.sv[0].replace('cnv', '')) + int(sample.sv[1].replace('cnv', ''))
		if total_cn < 4:
			return

	# allele fraction profile is not informative -- i.e. it's empty
	if sample.hap[0].af_mean_gene == -1:		
		sample.hap[0].add_dup(total_cn - 1)
		sample.gt = True
		return

	fit_maf, fit_cn = sample.hap[0].fit_data(total_cn, int(args_.gene_dict[args_.target_gene]['hg19_start']), int(args_.gene_dict[args_.target_gene]['hg19_end']))
	sample.hap[0].add_dup(fit_cn)
	sample.hap[1].add_dup(total_cn - fit_cn)
	sample.gt = True

def	cyp2a6_svcomb(sample):
	if sample.gt:
		return
	gt = []
	for sv in sample.sv:
		if sv == 'cnv0':
			gt.append(args_.star_dict['*4'])
		elif sv == 'cnv2':
			gt.append(new_dup(args_.star_dict['*1'], sv))
		elif sv == 'gc_e1e2':
			gt.append(args_.star_dict['*12'])
		elif sv == 'gc_e1e4':
			gt.append(args_.star_dict['*34'])			
		elif sv == 'dup2':
			gt.append(new_tandem(sv, ['*1', '*S1']))			
		elif sv == 'gc_e9':
			gt.append(new_tandem(sv, ['*1', '*S2']))
		elif sv == 'dup7':
			gt.append(new_tandem(sv, ['*1', '*S3']))
		elif sv == 'dup7x2':
			gt.append(new_tandem(sv, ['*1', '*S3', '*S3']))
		elif sv == 'dup7b':
			gt.append(new_tandem(sv, ['*1', '*S6']))
	if len(gt) == 2:
		sample.hap[0].cand = [gt[0]]
		sample.hap[1].cand = [gt[1]]
		sample.gt = True	
																		
def	cyp2d6_svcomb(sample):
	if sample.gt:
		return
	gt = []
	for sv in sample.sv:
		if sv == 'cnv0':
			gt.append(args_.star_dict['*5'])
		elif sv == 'gc_i1e9' and which_has(sample, ['*68', '*4']):
			gt.append(new_tandem(sv, ['*68', '*4']))
		elif sv == 'gc_i1e9' and which_has(sample, ['*S1', '*1']):
			gt.append(new_tandem(sv, ['*S1', '*1']))				
		elif sv == 'gc_e9' and which_has(sample, ['*4N', '*4']):
			gt.append(new_tandem(sv, ['*4N', '*4']))
		elif sv == 'gc_e9' and which_has(sample, ['*36', '*10']):
			gt.append(new_tandem(sv, ['*36', '*10']))
		elif sv == 'gc_e9' and which_has(sample, ['*83', '*2']):
			gt.append(new_tandem(sv, ['*83', '*2']))
		elif sv == 'gc_7to6_i4' and which_has(sample, ['13A', '*2']): 
			gt.append(new_tandem(sv, ['13A', '*2']))			
		elif sv == 'gc_7to6_i1':
			gt.append(args_.star_dict['*13B'])
		elif sv == 'gc_e1e7':
			gt.append(args_.star_dict['*13C'])
	cnv = None
	for sv in sample.sv:
		if 'cnv' in sv and sv != 'cnv0':
			cnv = sv
	if cnv:
		if '*68+*4' in [x.name for x in gt]:
			svcomb_tandem_cnv(gt, sample, ['*68', '*4'], cnv)
		elif '*S1+*1' in [x.name for x in gt]:
			svcomb_tandem_cnv(gt, sample, ['*S1', '*1'], cnv)			
		elif '*4N+*4' in [x.name for x in gt]:
			svcomb_tandem_cnv(gt, sample, ['*4N', '*4'], cnv)
		elif '*36+*10' in [x.name for x in gt]:
			svcomb_tandem_cnv(gt, sample, ['*36', '*10'], cnv)
		elif '*83+*2' in [x.name for x in gt]:
			svcomb_tandem_cnv(gt, sample, ['*83', '*2'], cnv)	
		elif '*13A+*2' in [x.name for x in gt]:
			svcomb_tandem_cnv(gt, sample, ['*13A', '*2'], cnv)
		elif '*13B' in [x.name for x in gt]:
			svcomb_sv1_cnv(gt, sample, '*13B', cnv)
	if len(gt) == 2:
		sample.hap[0].cand = [gt[0]]
		sample.hap[1].cand = [gt[1]]
		sample.gt = True

def svcomb_sv1_cnv(gt, sample, sX_name, cnv):
	i = which_has(sample, [sX_name])
	if not i:
		return
	if i != 3:
		j = {0: 1, 1: 0}[i - 1]
	elif i == 3:
		j = 0
	l = copy.deepcopy(sample.hap[j].cand)
	remove_sv(l)
	sY = new_dup(l[0], cnv)
	gt.insert(j, sY)
			
def svcomb_tandem_cnv(gt, sample, tandem, cnv):
	i = which_has(sample, [tandem[0], tandem[1]])
	if not i:
		return
	if i != 3:
		j = {0: 1, 1: 0}[i - 1]
		l = copy.deepcopy(sample.hap[j].cand)
		remove_sv(l)
		sX = new_dup(l[0], cnv)
		gt.insert(j, sX)
	elif i == 3:
		for x in sample.hap[0].cand:
			if x.name == tandem[1]:
				gt.append(new_dup(x, cnv))
				break
						
def	gstt1_svcomb(sample):
	if sample.gt:
		return
	gt = []
	for sv in sample.sv:
		if sv == 'cnv0':
			gt.append(args_.star_dict["*2"])
	if len(gt) == 2:
		sample.hap[0].cand = [gt[0]]
		sample.hap[1].cand = [gt[1]]		
		sample.gt = True

def	gstm1_svcomb(sample):
	if sample.gt:
		return
	gt = []
	for sv in sample.sv:
		if sv == 'cnv0':
			gt.append(args_.star_dict["*2"])
	if len(gt) == 2:
		sample.hap[0].cand = [gt[0]]
		sample.hap[1].cand = [gt[1]]		
		sample.gt = True	

def	ugt2b17_svcomb(sample):
	if sample.gt:
		return
	gt = []
	for sv in sample.sv:
		if sv == 'cnv0':
			gt.append(args_.star_dict["*2"])
	if len(gt) == 2:
		sample.hap[0].cand = [gt[0]]
		sample.hap[1].cand = [gt[1]]
		sample.gt = True		

##############################################################################

def remove_extra_s1():
	def f(l):
		if len(l) == 1:
			return
		for i in reversed(range(len(l))):
			if l[i].name == '*1':
				del l[i]

	for name, sample in args_.samples.items():
		f(sample.hap[0].cand)
		f(sample.hap[1].cand)
		f(sample.dip_cand)

def write_result_file():
	float2str = lambda x: '.' if x == -1 else '{:.2f}'.format(x)
	list2str = lambda x: '.' if not x else ','.join([str(x) for x in x])

	with open(args_.output_prefix + '.stargazer-genotype.txt', 'w') as f:
		header = ['name', 'status', 'hap1_main', 'hap2_main', 'hap1_cand', 'hap2_cand', 'hap1_score', 'hap2_score', 'dip_score', 'phenotype', 'dip_sv', 'hap1_sv', 'hap2_sv', 'ssr', 'dip_cand', 'hap1_main_core', 'hap2_main_core', 'hap1_main_tag', 'hap2_main_tag', 'hap1_af_mean_gene', 'hap2_af_mean_gene', 'hap1_af_mean_main', 'hap2_af_mean_main']
		f.write('\t'.join(header) + '\n')
		for name in sorted(args_.samples):
			sample = args_.samples[name]
			fields = ['.' for x in header]
			status = 'g' if sample.gt else 'ng'
			fields[header.index('name')] = name	
			fields[header.index('status')] = status
			if status == 'g': fields[header.index('hap1_main')] = sample.hap[0].cand[0].name
			if status == 'g': fields[header.index('hap2_main')] = sample.hap[1].cand[0].name
			if status != 'qc': fields[header.index('hap1_cand')] = ','.join([x.name for x in sample.hap[0].cand])
			if status != 'qc': fields[header.index('hap2_cand')] = ','.join([x.name for x in sample.hap[1].cand])
			if status == 'g': fields[header.index('hap1_score')] = str(sample.hap[0].cand[0].score)
			if status == 'g': fields[header.index('hap2_score')] = str(sample.hap[1].cand[0].score)
			if status == 'g': fields[header.index('dip_score')] = str(sample.dip_score)
			if status == 'g': fields[header.index('phenotype')] = sample.pt
			if status != 'qc': fields[header.index('dip_sv')] = ','.join(sample.sv)
			if status == 'g': fields[header.index('hap1_sv')] = sample.hap[0].sv
			if status == 'g': fields[header.index('hap2_sv')] = sample.hap[1].sv
			if status != 'qc': fields[header.index('ssr')] = sample.ssr
			if status != 'qc': fields[header.index('dip_cand')] = ','.join([x.name for x in sample.dip_cand])
			if status != 'qc': fields[header.index('hap1_main_core')] = list2str([x.summary() for x in sample.hap[0].obs if x in sample.hap[0].cand[0].core])
			if status != 'qc': fields[header.index('hap2_main_core')] = list2str([x.summary() for x in sample.hap[1].obs if x in sample.hap[1].cand[0].core])
			if status != 'qc': fields[header.index('hap1_main_tag')] = list2str([x.summary() for x in sample.hap[0].obs if x in sample.hap[0].cand[0].tag])
			if status != 'qc': fields[header.index('hap2_main_tag')] = list2str([x.summary() for x in sample.hap[1].obs if x in sample.hap[1].cand[0].tag])			
			if status != 'qc': fields[header.index('hap1_af_mean_gene')] = float2str(sample.hap[0].af_mean_gene(helper.get_hg19_start(args_.target_gene), helper.get_hg19_end(args_.target_gene)))
			if status != 'qc': fields[header.index('hap2_af_mean_gene')] = float2str(sample.hap[1].af_mean_gene(helper.get_hg19_start(args_.target_gene), helper.get_hg19_end(args_.target_gene)))
			if status != 'qc': fields[header.index('hap1_af_mean_main')] = float2str(sample.hap[0].af_mean_main)
			if status != 'qc': fields[header.index('hap2_af_mean_main')] = float2str(sample.hap[1].af_mean_main)
			f.write('\t'.join(fields) + '\n')

def predict_phenotypes():
	operators = {'<': operator.lt, '<=': operator.le, '>': operator.gt, '>=': operator.ge, '==': operator.eq}
	phenotypes = {}

	with open(args_.program_dir + '/phenotype_table.txt') as f:
		header = next(f).strip().split('\t')
		for line in f:
			fields = line.strip().split('\t')
			gene = fields[header.index('gene')]
			if gene != args_.target_gene:
				continue
			name = fields[header.index('name')]
			rules = fields[header.index('rules')].strip(',').split(',')
			phenotypes[name] = rules

	for name, sample in args_.samples.items():
		if sample.dip_score == 'unknown':
			sample.pt = 'unknown'
			continue
		for phenotype, rules in phenotypes.items():
			found = True
			for rule in rules:
				op = rule.split(':')[0]
				score = float(rule.split(':')[1])
				if not operators[op](sample.dip_score, score):
					found = False
					break
			if found:
				sample.pt = phenotype
				break

def order_haplotypes():
	for name, sample in args_.samples.items():
		if not sample.gt:
			continue
		if helper.sort_star_names([sample.hap[0].cand[0].name, sample.hap[1].cand[0].name])[0] == sample.hap[1].cand[0].name:
			sample.hap[0], sample.hap[1] = sample.hap[1], sample.hap[0]

def log_mode():
	with open(args_.log, 'a') as f:
		f.write(f'\n{args_.line_break}\n')
		f.write('Step 1/9: Determining genotype mode...\n\n')
		f.write('Status: Completed\n\n')		
		f.write(f'Target gene: {args_.target_gene.upper()}\n')
		f.write('Target paralog: {}\n'.format('N/A' if helper.get_paralog(args_.target_gene) == '.' else helper.get_paralog(args_.target_gene).upper()))
		f.write(f'Target region: chr{args_.target_region}\n\n')
		f.write('Control type: {}\n'.format('Control gene' if args_.control_gene else 'Copy number-stable region' if args_.control_type == 'cnsr' else 'Custom locus'))
		f.write('Control gene: {}\n'.format('N/A' if not args_.control_gene else args_.control_gene.upper()))
		f.write(f'Control region: chr{args_.control_region}\n\n')
		f.write('Input data source: {}\n'.format('Targeted sequencing' if args_.data_type == 'ts' else 'Whole genome sequencing' if args_.data_type == 'wgs' else 'Single nucleotide polymorphism array'))
		f.write(f'VCF-only mode is on: {args_.vcf_only}\n')
		f.write(f'Imputation mode is on: {args_.impute}\n\n')
		f.write(f'Enzyme function: {helper.get_function(args_.target_gene)}\n')
		f.write(f'PharmVar member: {helper.get_pv_member(args_.target_gene)}\n')
		f.write(f'DPSV member: {helper.get_dpsv_member(args_.target_gene)}\n')
		f.write(f'\n{args_.line_break}\n')

def assess_vcf(input_vcf):
	if not args_.vcf_only:
		# check whether the sample list is identical between VCF and GDF
		with open(args_.gdf) as f:
			gdf_samples = [x.replace('Depth_for_', '') for x in f.readline().strip().split('\t')[3:]]	
		vcf_samples = input_vcf.header[9:]	
		if len(gdf_samples) != len(vcf_samples):
			raise TypeError(f'The sample size differs between the VCF file (N={len(vcf_samples)}) and the GDF file (N={len(gdf_samples)})')
		if len(set(gdf_samples + vcf_samples)) != len(vcf_samples):
			raise TypeError(f'Two different sets of samples were detected from the VCF file and the GDF file')
		for i in range(len(vcf_samples)):
			if vcf_samples[i] != gdf_samples[i]:
				raise TypeError(f"The order of samples differs between the VCF file ('{vcf_samples[i]}') and the GDF file ('{gdf_samples[i]}') at sample index {i}")

		# make sure the sample size > 1 when using TS data
		if len(vcf_samples) < 5 and args_.data_type == 'ts':
			raise TypeError(f"Genotyping with TS data requires at least five samples (the current sample size is {len(vcf_samples)})")

	log_dict = {'row': 0, 'AD': 0, 'phased': 0, 'unphased': 0, 'both': 0}

	for fields in input_vcf.data:
		chrom = fields[0].replace('chr', '')
		pos = fields[1]
		format = fields[8].split(':')

		# Check GT field
		if 'GT' not in format:
			ValueError('GT field not found [{}]'.format(pos))

		# Check AD field
		if 'AD' in format:
			log_dict['AD'] += 1

		# Check phasing status
		def f(x):
			gt = x.split(':')[format.index('GT')]
			if '/' in gt:
				return '/'
			elif '|' in gt:
				return '|'
			else:
				if chrom == 'X' or chrom == 'Y':
					return
				else:
					raise ValueError('Genotype separator not found for autosomal chromosome chr{}:{} GT=[{}]'.format(chrom, pos, gt))

		separators = set([f(x) for x in fields[9:] if f(x)])

		log_dict['row'] += 1
		if len(separators) == 1:
			if '|' in separators:
				log_dict['phased'] += 1
			else:
				log_dict['unphased'] += 1
		else:
			log_dict['both'] += 1

	# Check if input VCF is empty
	if log_dict['row'] == 0:
		args_.vcf_empty = True

	# Determine AD mode
	if log_dict['row'] > 0 and log_dict['AD'] / log_dict['row'] > 0.8:
		args_.vcf_ad = True

	# Determine phasing mode
	if log_dict['phased'] == log_dict['row']:
		args_.vcf_sep = '|'
	elif log_dict['unphased'] == log_dict['row']:
		args_.vcf_sep = '/'
	else:
		args_.vcf_sep = 'b'

	with open(args_.log, 'a') as f:	
		f.write('Step 2/9: Assessing input VCF...\n\n')
		f.write('Status: Completed\n\n')
		f.write('Samples total: {}\n'.format(len(input_vcf.header[9:])))
		f.write('Markers total: {}\n\n'.format(log_dict['row']))
		f.write('Markers with allelic depth: {}\n\n'.format(log_dict['AD']))
		f.write('Markers unphased: {}\n'.format(log_dict['unphased']))
		f.write('Markers phased: {}\n'.format(log_dict['phased']))
		f.write('Markers partially phased: {}\n'.format(log_dict['both']))
		f.write(f'\n{args_.line_break}\n')

def process_vcf(input_vcf):
	log_dict = {'IA': 0, 'allelic_imbalance': 0, 's50': 0}
	processed_vcf = helper.copy_vcf(input_vcf, ['header'])	
	processed_vcf.meta = [
		'##fileformat=VCFv4.2\n',
		'##fileDate={}\n'.format(args_.date),
		'##source={}\n'.format(args_.stargazer_version),
		'##reference=hg19\n',
		'##INFO=<ID=FE,Number=A,Type=String,Description="Functional Effect">\n',
		'##INFO=<ID=PS,Number=1,Type=String,Description="Phasing Status (A, in preparation; B1, ready for phasing as is; B2, ready for phasing after conformation to reference VCF; C1, excluded from phasing because marker is absent in reference VCF; C2, excluded from phasing because marker has different REF allele; C3, excluded from phasing because marker has no overlapping ALT alleles; D1, statistically phased; D2, manually phased with certainty; D3, manually phased without certainty; D4, already phased; D5, manually phased by extension; E, omitted during statistical phasing)">\n',
		'##INFO=<ID=RV,Number=A,Type=String,Description="Reverting Variation">\n',
		'##INFO=<ID=SO,Number=A,Type=String,Description="=Sequence Ontology">\n',
		'##INFO=<ID=VI,Number=A,Type=String,Description="Variant Impact">\n',
		'##FILTER=<ID=IA,Description="Invalid Allele">\n',
		'##FILTER=<ID=s50,Description="Less than 50% of samples have data">\n',
		'##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n',
		'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
		'##FORMAT=<ID=HE,Number=4,Type=Integer,Description="Matching scores computed by the phase-by-extension algorithm">\n'
	]

	for fields in input_vcf.data:
		chr, ref, alt, fmt, flt, inf = fields[0].replace('chr', ''), fields[3], fields[4].split(','), fields[8].split(':'), [], ['PS=D4'] if args_.vcf_sep == '|' else ['PS=A']

		def f(x):
			gt_field = x.split(':')[fmt.index('GT')]
			
			# Determine genotype separator
			if '/' in gt_field:
				gt_sep = '/'
			elif '|' in gt_field:
				gt_sep = '|'
			else:
				gt_sep = ''

			# Unphase genotype if input VCF is partially phased
			if args_.vcf_sep == 'b' and gt_sep == '|':
				if gt_field == '.|.':
					gt_field = './.'
				else:
					gt_field = '/'.join(sorted(gt_field.split('|'), key = lambda x: int(x)))

			# Conform genotype for sex chromosomes if necessary
			if not gt_sep and (chr == 'X' or chr == 'Y'):
				gt_field = '0|' + gt_field if args_.vcf_sep == '|' else '0/' + gt_field

			# Get AD field information
			if args_.vcf_ad:
				ad_field = x.split(':')[fmt.index('AD')]
				
				if gt_field == './.':
					ad_field = ','.join(['0'] * (len(alt) + 1))
				
				# Some tools such as Genalice produce variable-length AD fields (i.e., single AD value if sample is homozygous for ALT allele)
				if len(ad_field.split(',')) == 1 and gt_sep and gt_field.split(gt_sep)[0] == gt_field.split(gt_sep)[1] and gt_field.split(gt_sep)[0] != '0' and len(alt) == 1:
					if 'DP' in fmt:
						dp_field = x.split(':')[fmt.index('DP')]
						ad_field = str(int(dp_field) - int(ad_field)) + ',' + ad_field
					else:
						ad_field = '0,' + ad_field
				ad_field = ':' + ad_field

			else:
				ad_field = ''

			return gt_field + ad_field

		fields[9:] = [f(x) for x in fields[9:]]

		# Check invalid allele
		if ref == 'I' or ref == 'D' or '.' in alt or 'D' in alt:
			flt.append('IA')
			log_dict['IA'] += 1

		# Check high missingness	
		if ['.' in x.split(':')[0] for x in fields[9:]].count(True) / len(fields[9:]) > 0.5:
			flt.append('s50')
			log_dict['s50'] += 1

		# Define quick function 3
		def qf3(x):
			gt = x.split(':')[0]
			if '.' in gt:
				return False
			if '|' in gt:
				gt = gt.split('|')
			else:
				gt = gt.split('/')
			if gt[0] == gt[1]:
				return False				
			ad = [int(y) for y in x.split(':')[1].split(',')]
			if sum(ad) == 0:
				return False
			return max(ad) / sum(ad)

		# Check allelic imbalance
		if args_.vcf_ad:
			ratios = [qf3(x) for x in fields[9:] if qf3(x)]
			if ratios:
				median = statistics.median(ratios)
				if median > 0.8 or median < 0.2:
					log_dict['allelic_imbalance'] += 1				

		fields[5] = '.'
		fields[6] = ';'.join(flt) if flt else 'PASS'
		fields[7] = ';'.join(inf)
		fields[8] = 'GT:AD' if args_.vcf_ad else 'GT'
		processed_vcf.data.append(fields)

	with open(args_.log, 'a') as f:
		f.write('Step 3/9: Processing input VCF...\n\n')
		f.write('Status: Completed\n\n')
		f.write('Markers with allelic imbalance: {}\n'.format(log_dict['allelic_imbalance']))
		f.write('Markers with high missingness: {}\n'.format(log_dict['s50']))
		f.write('Markers with invalid allele: {}\n'.format(log_dict['IA']))
		f.write(f'\n{args_.line_break}\n')

	return processed_vcf

def conform_vcf(processed_vcf, ref_vcf):
	log_dict = {'status': 'Skipped because input VCF is empty ', 'row': 0, 'filtered': 0, 'PS=B1': 0, 'PS=B2': 0, 'PS=C1': 0, 'PS=C2': 0, 'PS=C3': 0}

	def write_log():
		with open(args_.log, 'a') as f:
			f.write('Step 4/14: Conforming input VCF...\n\n')
			f.write('Status: {}\n\n'.format(log_dict['status']))
			f.write('Markers total: {}\n'.format(len(processed_vcf.data)))
			f.write('Markers filtered: {}\n'.format(log_dict['filtered']))
			f.write('Markers remaining: {}\n\n'.format(log_dict['row']))
			f.write('Markers phasable: {}\n'.format(log_dict['PS=B1'] + log_dict['PS=B2']))
			f.write('Markers ready: {}\n'.format(log_dict['PS=B1']))
			f.write('Markers conformed: {}\n\n'.format(log_dict['PS=B2']))
			f.write('Markers unphasable: {}\n'.format(log_dict['PS=C1'] + log_dict['PS=C2'] + log_dict['PS=C3']))
			f.write('Markers absent in reference VCF: {}\n'.format(log_dict['PS=C1']))
			f.write('Markers with different REF allele: {}\n'.format(log_dict['PS=C2']))
			f.write('Markers with no overlapping ALT alleles: {}\n'.format(log_dict['PS=C3']))
			f.write(f'\n{args_.line_break}\n')

	if args_.vcf_empty:
		write_log()
		return processed_vcf

	if args_.vcf_sep == '|':
		log_dict['status'] = 'Skipped because input VCF is already fully phased'
		write_log()
		return processed_vcf
		
	log_dict['status'] = 'Completed'
	conformed_vcf = helper.copy_vcf(processed_vcf, ['meta', 'header'])

	for fields1 in processed_vcf.data:
		chr1, pos1, ref1, alt1, flt1, inf1 = fields1[0], fields1[1], fields1[3], fields1[4].split(','), fields1[6], fields1[7].split(';')

		if flt1 != 'PASS':
			log_dict['filtered'] += 1
			continue

		log_dict['row'] += 1	
	
		is_found = False
	
		for i in range(len(ref_vcf.data)):
			fields2 = ref_vcf.data[i]
			chr2, pos2, ref2, alt2 = fields2[0], fields2[1], fields2[3], fields2[4].split(',')
	
			# Keep looking if not found
			if chr1 != chr2 or pos1 != pos2:
				continue

			# Although found, will not be phased
			if ref1 != ref2:
				fields3 = ref_vcf.data[i + 1]
				pos3 = fields3[1]
				ref3 = fields3[3]
			
				# Check if the next line matches
				if pos1 == pos3 and ref1 == ref3:
					continue
				else:
					log_dict['PS=C2'] += 1
					inf1 = ['PS=C2' if x == 'PS=A' else x for x in inf1]
					is_found = True
					break
		
			# There are no overlapping ALT alleles
			if not list(set(alt1) & set(alt2)):
				is_found = True
				log_dict['PS=C3'] += 1
				inf1 = ['PS=C3' if x == 'PS=A' else x for x in inf1]
				break

			# Found and perfectly matched, no need to conform
			if alt1 == alt2 or len(set(alt1) & set(alt2)) == 0:
				is_found = True
				log_dict['PS=B1'] += 1
				inf1 = ['PS=B1' if x == 'PS=A' else x for x in inf1]
				break

			# Although found, missing one or more ALT alleles
			if set(alt1).issubset(alt2) and len(alt2) > len(alt1):
				diff = len(alt2) - len(alt1)
				for allele in alt2:
					if allele not in alt1:
						alt1.append(allele)
				if args_.vcf_ad:
					fields1[9:] = [x + ',0' * diff for x in fields1[9:]]

			# Although same ALT alleles, wrong order
			if set(alt1) == set(alt2):
				is_found = True
				mapping = {0: 0}
				for i in range(len(alt1)):
					mapping[i + 1] = alt2.index(alt1[i]) + 1
				fields1[4] = ','.join(alt2) # update ALT alleles
				log_dict['PS=B2'] += 1
				inf1 = ['PS=B2' if x == 'PS=A' else x for x in inf1]
				
				def f(x):
					gt = x.split(':')[0].split('/')
					for i in [0, 1]:
						if gt[i] != '0' and gt[i] != '.':
							gt[i] = str(mapping[int(gt[i])])
					if not args_.vcf_ad:
						return '/'.join(gt)				
					ad1 = x.split(':')[1].split(',')
					ad2 = [0 for y in ad1]
					for i in range(len(ad2)):
						ad2[mapping[i]] = ad1[i]
					return '/'.join(gt) + ':' + ','.join(ad2)				

				fields1[9:] = [f(x) for x in fields1[9:]]
				break

		if not is_found:
			inf1 = ['PS=C1' if x == 'PS=A' else x for x in inf1]
			log_dict['PS=C1'] += 1	
	
		fields1[7] = ';'.join(inf1)
		conformed_vcf.data.append(fields1)

	write_log()
	return conformed_vcf

def phase_vcf(conformed_vcf):
	log_dict = {'status': 'Skipped because input VCF is empty', 'attempted': 0, 'PS=D1': 0, 'PS=E': 0}

	def write_log():
		with open(args_.log, 'a') as f:
			f.write('Step 5/9: Statistically phasing input VCF...\n\n')
			f.write('Status: {}\n\n'.format(log_dict['status']))
			f.write('Markers attempted: {}\n'.format(log_dict['attempted']))
			f.write('Markers phased: {}\n'.format(log_dict['PS=D1']))
			f.write('Markers omitted: {}\n'.format(log_dict['PS=E']))
			f.write(f'\n{args_.line_break}\n')

	if args_.vcf_empty:
		write_log()
		return conformed_vcf

	if args_.vcf_sep == '|':
		log_dict['status'] = 'Skipped because input VCF is already fully phased'
		write_log()
		return conformed_vcf
	
	phaseme_vcf = helper.copy_vcf(conformed_vcf, ['meta', 'header', 'data'])

	for i in reversed(range(len(phaseme_vcf.data))):
		fields = phaseme_vcf.data[i]
		inf = fields[7].split(';')
		if any(['PS=C' in x for x in inf]):
			del phaseme_vcf.data[i]

	log_dict['attempted'] = len(phaseme_vcf.data)

	if log_dict['attempted'] == 0:
		log_dict['status'] = 'Skipped because none of the markers are eligible for phasing'
		write_log()
		return conformed_vcf

	if log_dict['attempted'] == 1:
		log_dict['status'] = 'Skipped because there is only single marker eligible for phasing (manually phased)'
		combined_vcf = helper.copy_vcf(conformed_vcf, ['meta', 'header', 'data'])
		
		for i in range(len(combined_vcf.data)):
			fields = combined_vcf.data[i]
			inf = fields[7].split(';')
			if not any(['PS=C' in x for x in inf]):
				fields[7] = 'PS=D2'
				combined_vcf.data[i] = fields[:9] + [x.replace('/', '|') for x in fields[9:]]
				break
		write_log()
		return combined_vcf

	log_dict['status'] = 'Completed'
	helper.write_vcf(phaseme_vcf, f'{args_.extended_prefix}.project/phaseme.vcf')
	subprocess.call(['java', '-Xmx2g', '-jar', args_.program_dir + '/beagle.25Nov19.28d.jar', f'gt={args_.extended_prefix}.project/phaseme.vcf', 'chrom=' + args_.target_region, 'ref=' + args_.ref_vcf, f'out={args_.extended_prefix}.project/phased', 'impute=' + str(args_.impute).lower()], stdout = open(f'{args_.extended_prefix}.project/phased.log', 'a'), stderr = open(f'{args_.extended_prefix}.project/phased.log', 'a'))
	subprocess.call(['gunzip', f'{args_.extended_prefix}.project/phased.vcf.gz'])
	phased_vcf = helper.read_vcf_simple(f'{args_.extended_prefix}.project/phased.vcf')	
	combined_vcf = helper.copy_vcf(conformed_vcf, ['meta', 'header'])
	
	for fields1 in conformed_vcf.data:
		pos1, ref1, alt1, inf1 = fields1[1], fields1[3], fields1[4], fields1[7].split(';')

		if any(['PS=C' in x for x in inf1]):
			fields1[9:] = [x.replace('./.', '0/0') for x in fields1[9:]]

			if all([x.split(':')[0].split('/')[0] == x.split(':')[0].split('/')[1] for x in fields1[9:]]):
				fields1[9:] = [x.replace('/', '|') for x in fields1[9:]]
				fields1[7] = 'PS=D2'

			combined_vcf.data.append(fields1)
			continue

		is_found = False

		for fields2 in phased_vcf.data:			
			pos2, ref2, alt2 = fields2[1], fields2[3], fields2[4]

			gt = lambda x: fields2[9:][x].split(':')[0]
			ad = lambda x: fields1[9:][x].split(':')[1]
			idx = list(range(len(fields1[9:])))

			if pos1 == pos2 and ref1 == ref2 and alt1 == alt2:
				is_found = True
				fields1[9:] = [gt(i) + ':' + ad(i) for i in idx] if args_.vcf_ad else [gt(i) for i in idx]
				log_dict['PS=D1'] += 1
				inf1 = ['PS=D1' if 'PS=B' in x else x for x in inf1]
				break

		if not is_found:
			log_dict['PS=E'] += 1
			inf1 = ['PS=E' if 'PS=B' in x else x for x in inf1]

		fields1[7] = ';'.join(inf1)
		combined_vcf.data.append(fields1)

	write_log()
	return combined_vcf

def annotate_vcf(combined_vcf):
	log_dict = {'stargazer_membership': 0, 'low_impact': 0, 'moderate_impact': 0, 'high_impact': 0, 'reverting_variation': 0}
	annotated_vcf = helper.copy_vcf(combined_vcf, ['meta', 'header', 'data'])
	undetected_revertants = [x for x in args_.snp_list if x.rev]
	
	for fields in annotated_vcf.data:
		pos, ref, alt, inf = fields[1], fields[3], fields[4].split(','), fields[7].split(';')
		vi_list = []; rv_list = []; so_list = []; fe_list = []

		for var in alt:
			filtered = [x for x in args_.snp_list if pos == x.pos and ref == x.hg and (var == x.var or var == x.wt)]
			if filtered:
				vi, rv, so, fe = filtered[0].impact, 'revertant_true' if filtered[0].rev else 'revertant_false', filtered[0].so, filtered[0].effect
				log_dict[vi] += 1
				undetected_revertants = [x for x in undetected_revertants if x != filtered[0]]
				if filtered[0].rev:
					log_dict['reverting_variation'] += 1
			else:
				vi, rv, so, fe = 'unknown_impact', 'revertant_unknown', 'unknown_variant', 'unknown_effect'
			vi_list.append(vi); rv_list.append(rv); so_list.append(so); fe_list.append(fe)

		if any([x != 'unknown_impact' for x in vi_list]):
			log_dict['stargazer_membership'] += 1

		inf.append('VI=' + ','.join(vi_list)); inf.append('RV=' + ','.join(rv_list)); inf.append('SO=' + ','.join(so_list)); inf.append('FE=' + ','.join(fe_list))
		fields[7] = ';'.join(inf)
	
	# Manually add variants that are part of the genome assembly
	if args_.data_type != 'chip':
		chr = helper.get_chr(args_.target_gene)
		dat, fmt = (['0|0:0,0' for x in annotated_vcf.header[9:]], 'GT:AD') if args_.vcf_ad else (['0|0' for x in annotated_vcf.header[9:]], 'GT')
		for snp in undetected_revertants:
			inf = ';'.join(['PS=D2', 'VI=' + snp.impact, 'RV=revertant_true', 'SO=' + snp.so, 'FE=' + snp.effect])
			fields = [chr, snp.pos, snp.id, snp.hg, snp.wt, '.', 'PASS', inf, fmt] + dat
			annotated_vcf.data.append(fields)
		annotated_vcf.data.sort(key = lambda x: int(x[1]))

	# Manually phase, without certainty, any unphased revertants
	for fields in annotated_vcf.data:
		inf = fields[7].split(';')
		if any(['PS=C' in x for x in inf]) and 'RV=revertant_true' in inf:
			fields[7] = ';'.join(['PS=D3' if 'PS=C' in x else x for x in inf])
			fields[9:] = [x.replace('/', '|') for x in fields[9:]]
	
	with open(args_.log, 'a') as f:
		f.write('Step 6/9: Annotating input VCF...\n\n')
		f.write('Status: Completed\n\n')		
		f.write('Markers with Stargazer membership: {}\n\n'.format(log_dict['stargazer_membership']))
		f.write('Variants with low impact: {}\n'.format(log_dict['low_impact']))
		f.write('Variants with moderate impact: {}\n'.format(log_dict['moderate_impact']))
		f.write('Variants with high impact: {}\n\n'.format(log_dict['high_impact']))
		f.write('Variants reverted to wild type: {}/{}\n'.format(log_dict['reverting_variation'], len([x for x in args_.snp_list if x.rev])))
		f.write(f'\n{args_.line_break}\n')

	return annotated_vcf

def account_vcf(annotated_vcf):
	accounted_vcf = helper.copy_vcf(annotated_vcf, ['meta', 'header', 'data'])
	accounted_vcf.meta = ['##reference=hg19-{}*1\n'.format(args_.target_gene.upper()) if '##reference=' in x else x for x in accounted_vcf.meta]

	for fields in accounted_vcf.data:
		ref, alt, inf = fields[3], fields[4].split(','), fields[7].split(';')
		rv_list = [x for x in inf if 'RV=' in x][0].replace('RV=', '').split(',')
		if 'revertant_true' not in rv_list:
			continue
		i = rv_list.index('revertant_true')
		fields[3] = alt[i]
		fields[4] = ','.join([ref if x == alt[i] else x for x in alt])

		def f(x):
			gt = x.split(':')[0].split('|')
			field = '|'.join([str(i + 1) if y == '0' else '0' if y == str(i + 1) else y for y in gt])
			if args_.vcf_ad:				
				ad = x.split(':')[1].split(',')
				field += ':' + ','.join([ad[i + 1] if y == 0 else ad[0] if y == i + 1 else ad[y] for y in range(len(ad))])
			return field	
		
		fields[9:] = [f(x) for x in fields[9:]]

	return accounted_vcf

def extend_vcf(accounted_vcf):
	log_dict = {'attempted': 0, 'phased': 0, 'omitted': 0}
	finalized_vcf = helper.copy_vcf(accounted_vcf, ['meta', 'header', 'data'])
	pseudo_samples = helper.vcf2samples(finalized_vcf)

	for fields in finalized_vcf.data:
	
		pos, ref, alt, inf, fmt = fields[1], fields[3], fields[4].split(','), fields[7].split(';'), fields[8]
		ps, vi = inf[0], inf[1]
		if 'PS=D' in ps or ('high_impact' not in vi and 'moderate_impact' not in vi):
			continue
		log_dict['attempted'] += 1

		def pbe(i):
			x = fields[i]
			gt = x.split(':')[0].split('/')
			if gt[0] == gt[1]:
				return x.replace('/', '|')			
			scores = [[0, 0], [0, 0]]
			for j in [0, 1]:
				if gt[j] == '0':
					continue
				idx = int(gt[j])
				target_snp = helper.SNP()
				target_snp.pos = pos
				target_snp.wt = ref
				target_snp.var = alt[idx - 1]
				relevant_stars = [v for k, v in args_.star_dict.items() if target_snp in v.core]
				name = finalized_vcf.header[i]
				for k in [0, 1]:
					for star in relevant_stars:
						score = 0
						for snp in pseudo_samples[name].hap[k].obs:
							if snp in star.core + star.tag:
								score += 1
						if score > scores[j][k]:
							scores[j][k] = score			
			a = scores[0]
			b = scores[1]
			flip = False
			if max(a) == max(b):
				if a[0] < a[1] and b[0] > b[1]:
					flip = True
				elif a[0] == a[1] and b[0] > b[1]:
					flip = True
				elif a[0] < a[1] and b[0] == b[1]:
					flip = True					
				else:
					pass
			else:
				if max(a) > max(b):
					if a[0] > a[1]:
						pass
					else:
						flip = True
				else:
					if b[0] > b[1]:
						flip = True
					else:
						pass
			if flip:
				result = f'{gt[1]}|{gt[0]}'
			else:
				result = f'{gt[0]}|{gt[1]}'
			if 'AD' in fmt:
				result = result + ':' + x.split(':')[1]
			result = result + f':{",".join([str(x) for x in a + b])}'
			return result

		new_fields = [pbe(i) for i in range(9, len(fields))]

		if not all(new_fields):
			log_dict['omitted'] += 1
			continue

		log_dict['phased'] += 1

		inf[0] = 'PS=D5'
		fields[7] = ';'.join(inf)
		fields[8] += ':HE'
		fields[9:] = new_fields

	with open(args_.log, 'a') as f:
		f.write('Step 7/9: Phasing input VCF by haplotype extension...\n\n')
		f.write('Status: Completed\n\n')
		f.write(f'Markers attempted: {log_dict["attempted"]}\n')
		f.write(f'Markers phased: {log_dict["phased"]}\n')
		f.write(f'Markers omitted: {log_dict["omitted"]}\n')
		f.write(f'\n{args_.line_break}\n')

	return finalized_vcf
				
def find_candidate_stars():
	for name, sample in args_.samples.items():
		f = lambda x: sorted([v for k, v in args_.star_dict.items() if set(v.core).issubset(x) and not (v.sv and v.sv not in sample.sv)], key = lambda x: x.rank)
		hap1_snp = [x for x in sample.hap[0].obs if x.wt != x.var]
		hap2_snp = [x for x in sample.hap[1].obs if x.wt != x.var]
		sample.hap[0].cand = f(hap1_snp)
		sample.hap[1].cand = f(hap2_snp)
		sample.dip_cand = f(list(set(hap1_snp + hap2_snp)))

def transfer_sv_data():
	if args_.vcf_only:
		for name, sample in args_.samples.items():
			sample.sv = ['no_sv', 'no_sv']
			sample.ssr = '.'
	else:
		for name, sample in args_.samples.items():
			with open(f'{args_.extended_prefix}.project/ssr/{name}.txt') as f:
				header = next(f).strip().split('\t')
				fields = next(f).strip().split('\t')
				sample.sv = [fields[header.index('seq1')], fields[header.index('seq2')]]
				sample.ssr = fields[header.index('ssr')]

def make_sv_calls():
	log_dict = {'status': 'Skipped because VCF-only mode is on'}
	def write_log():
		with open(args_.log, 'a') as f:
			f.write('Step 8/9: Detecting structural variants...\n\n')
			f.write('Status: {}\n'.format(log_dict['status']))
			f.write(f'\n{args_.line_break}\n')
	if args_.vcf_only:
		write_log()
		return
	log_dict['status'] = 'Completed'
	subprocess.call(['mkdir', f'{args_.extended_prefix}.project/ssr'])
	exit_code = subprocess.call(['Rscript', args_.program_dir + '/sv.R', args_.program_dir, args_.extended_prefix, args_.data_type, args_.target_gene, args_.control_gene, args_.control_region, args_.gdf, args_.control_type, args_.sample_list])
	if exit_code == 0:
		write_log()
	else:
		raise TypeError('Something bad happended during SV detection!')

def make_genotype_calls():
	for name, sample in args_.samples.items():
		if sample.sv == ['no_sv', 'no_sv']: sample.gt = True
		call_sv1(sample, 'cnv0', args_.deletion)
		call_cnv3(sample)
		call_cnv_plus(sample)
		callers = {'cyp2a6': cyp2a6, 'cyp2b6': cyp2b6, 'cyp2d6': cyp2d6, 'cyp2e1': cyp2e1, 'gstm1': gstm1, 'gstt1': gstt1, 'slc22a2': slc22a2, 'slco1b1': slco1b1, 'ugt1a4': ugt1a4, 'ugt2b15': ugt2b15, 'ugt2b17': ugt2b17}
		if args_.target_gene in callers:	
			callers[args_.target_gene](sample)

def plot_profiles():
	log_dict = {'status': 'Skipped because VCF-only mode is on'}

	def write_log():
		with open(args_.log, 'a') as f:
			f.write('Step 9/9: Plotting various profiles...\n')
			f.write('\n')
			f.write('Status: {}\n'.format(log_dict['status']))
			f.write(f'\n{args_.line_break}')

	if args_.vcf_only:
		write_log()
		return
		
	log_dict = {'status': 'Completed'}
	
	subprocess.call(['mkdir', f'{args_.extended_prefix}.project/af'])
	for name, sample in args_.samples.items():
		with open(f'{args_.extended_prefix}.project/af/{name}.txt', 'w') as f:
			header = ['pos', 'hap1_al', 'hap2_al', 'hap1_af', 'hap2_af']
			f.write('\t'.join(header) + '\n')
			filtered = [x for x in sample.hap[0].obs if x.td > 10 and x.het]
			if not filtered:
				continue
			for snp in filtered:
				hap1_snp = [x for x in sample.hap[0].obs if x.pos == snp.pos][0]
				hap2_snp = [x for x in sample.hap[1].obs if x.pos == snp.pos][0]
				fields = [snp.pos, hap1_snp.var, hap2_snp.var, str(hap1_snp.af), str(hap2_snp.af)]
				f.write('\t'.join(fields) + '\n')

	subprocess.call(['mkdir', f'{args_.extended_prefix}.project/plot'])
	exit_code = subprocess.call(['Rscript', args_.program_dir + '/plot.R', args_.extended_prefix, str(args_.detail)])	
	if exit_code == 0:
		write_log()
	else:
		raise TypeError('Something bad happended while plotting profiles!')	
	
def run(args):
	global args_
	args_ = args

	# assess the arguments
	if args_.gdf and args_.control_type == 'known' and not args_.control_gene: raise TypeError(f"Argument 'control_gene' not found, required because of arguments 'gdf' [={args_.gdf}] and 'control_type' [={args_.control_type}]")
	if args_.control_gene and not args_.gdf: raise TypeError(f"Argument 'gdf' not found, required because of argument 'control_gene' [={args_.control_gene}]")
	if args_.control_type == 'custom' and not args_.region: raise TypeError(f"Argument 'region' not found, required because of argument 'control_type' [={args_.control_type}]")
	if args_.control_type == 'cnsr' and helper.get_cnsr(args_.target_gene) == '.': raise TypeError("CNSR not avaialble as control for selected target gene, use either known control gene or custom locus")

	# update the container for the "genotype" tool
	args_.vcf_ad = False
	args_.vcf_sep = None
	args_.vcf_empty = False
	args_.vcf_only = (args_.control_type == 'known' and not args_.control_gene) or not args_.gdf
	args_.deletion = [v.name for k, v in args_.star_dict.items() if v.sv == 'cnv0'][0]
	args_.ref_vcf = '{}/1kgp_vcf/{}.vcf.gz'.format(args_.program_dir, args_.target_gene) if not args_.ref_vcf else args_.ref_vcf
	args_.control_gene = '.' if args_.control_type != 'known' else args_.control_gene
	args_.sample_list = ','.join(args_.sample_list) if args_.sample_list else '.'
	
	helper.create_dir(f'{args_.extended_prefix}.project')
		
	log_mode() # step 1
	ref_vcf = helper.read_vcf_simple(args_.ref_vcf)
	input_vcf = helper.read_vcf_region(args_.vcf, args_.target_region)
	assess_vcf(input_vcf) # step 2

	processed_vcf = process_vcf(input_vcf) # step 3
# 	helper.write_vcf(processed_vcf, f'{args_.extended_prefix}.project/processed.vcf')
	conformed_vcf = conform_vcf(processed_vcf, ref_vcf) # step 4
# 	helper.write_vcf(conformed_vcf, f'{args_.extended_prefix}.project/conformed.vcf')
	combined_vcf = phase_vcf(conformed_vcf) # step 5
# 	helper.write_vcf(combined_vcf, f'{args_.extended_prefix}.project/combined.vcf')
	annotated_vcf = annotate_vcf(combined_vcf) # step 6
# 	helper.write_vcf(annotated_vcf, f'{args_.extended_prefix}.project/annotated.vcf')
	accounted_vcf = account_vcf(annotated_vcf)
# 	helper.write_vcf(accounted_vcf, f'{args_.extended_prefix}.project/accounted.vcf')
	finalized_vcf = extend_vcf(accounted_vcf) # step 7
	helper.write_vcf(finalized_vcf, f'{args_.extended_prefix}.project/finalized.vcf')

	args_.samples = helper.vcf2samples(finalized_vcf)

	make_sv_calls() # step 8
	transfer_sv_data()

	find_candidate_stars()
	make_genotype_calls()
	remove_extra_s1()
	order_haplotypes()
	predict_phenotypes()
	write_result_file()
	plot_profiles() # step 9
