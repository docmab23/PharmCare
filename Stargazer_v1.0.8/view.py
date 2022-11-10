import os
import subprocess
import helper

def summary():
	alleles = {}
	allele_scores = {}
	hap_scores = {}
	dip_scores = {}	
	phenotypes = {}
	sv_counts = {}
	sv_types = {}
	s_sv = 0
	with open(args_.genotype) as f:
		header = next(f).strip().split('\t')
		for line in f:
			fields = line.strip().split('\t')
			hap1_main = fields[header.index('hap1_main')]
			hap2_main = fields[header.index('hap2_main')]
			hap1_score = fields[header.index('hap1_score')]
			hap2_score = fields[header.index('hap2_score')]
			hap1_sv = fields[header.index('hap1_sv')]
			hap2_sv = fields[header.index('hap2_sv')]
			dip_score = fields[header.index('dip_score')]
			phenotype = fields[header.index('phenotype')]
			dip_sv = fields[header.index('dip_sv')]

			if hap1_main not in alleles:
				alleles[hap1_main] = [hap1_sv, hap1_score, 0]
			if hap2_main not in alleles:
				alleles[hap2_main] = [hap2_sv, hap2_score, 0]
			if dip_score not in dip_scores:
				dip_scores[dip_score]= 0
			if phenotype not in phenotypes:
				phenotypes[phenotype] = 0
			if dip_sv != 'no_sv,no_sv':
				s_sv += 1
			if 'unknown' in dip_sv:
				raise ValueError
			else:
				sv_count = 2 - dip_sv.split(',').count('no_sv')	
			if sv_count not in sv_counts:
				sv_counts[sv_count] = 0
			if hap1_score not in hap_scores:
				hap_scores[hap1_score] = 0
			if hap2_score not in hap_scores:
				hap_scores[hap2_score] = 0
			if hap1_sv not in sv_types:
				sv_types[hap1_sv] = 0
			if hap2_sv not in sv_types:
				sv_types[hap2_sv] = 0				
			hap_scores[hap1_score] += 1
			hap_scores[hap2_score] += 1
			alleles[hap1_main][2] += 1
			alleles[hap2_main][2] += 1
			dip_scores[dip_score] += 1
			phenotypes[phenotype] += 1	
			sv_counts[sv_count] += 1
			sv_types[hap1_sv] += 1
			sv_types[hap2_sv] += 1
	a_total = sum([x[2] for x in alleles.values()])	
	s_total = a_total / 2
	
	# Count allele scores
	for allele in alleles:
		score = alleles[allele][1]
		if score not in allele_scores:
			allele_scores[score] = 0
		allele_scores[score] += 1
		
	header = ['type', 'name', 'sv', 'hap_score', 'count', 'percentage']
	
	
	with open(f'{args_.extended_prefix}.txt', 'w') as f:		
		
		def w(l):
			if len(l) != len(header):
				raise ValueError()
			f.write('\t'.join([str(x) for x in l]) + '\n')

		w(header)
		w(['samples', 'total', '.', '.', s_total, s_total/s_total])
		w(['samples', 'sv', '.', '.', s_sv, s_sv/s_total])
		for sv_count in sorted(sv_counts):
			w(['sv_count', sv_count, '.', '.', sv_counts[sv_count], sv_counts[sv_count] / s_total])
		w(['haps', 'total', '.', '.', a_total, a_total/a_total])		
		w(['haps', 'unique', '.', '.', len(alleles), '.'])
		for sv_type in sorted(sv_types):
			w(['sv_types', sv_type, '.', '.', sv_types[sv_type], sv_types[sv_type] / a_total])
		for allele in helper.sort_star_names(list(alleles)):
			sv = alleles[allele][0]
			hap_score = alleles[allele][1]
			n = alleles[allele][2]
			p = n / a_total			
			w(['alleles', allele, sv, hap_score, n, p])
		for allele_score in sorted(allele_scores):
			n = allele_scores[allele_score]
			p = n / len(alleles)
			w(['allele_scores', allele_score, '.', '.', n, p])
		for hap_score in sorted(hap_scores):
			n = hap_scores[hap_score]
			p = n / a_total
			w(['hap_scores', hap_score, '.', '.', n, p])
		for dip_score in sorted(dip_scores):
			n = dip_scores[dip_score]
			p = n / s_total
			w(['dip_scores', dip_score, '.', '.', n, p])
		for phenotype in sorted(phenotypes, key = helper.PHENOTYPES.index):
			n = phenotypes[phenotype]
			p = phenotypes[phenotype] / s_total
			w(['phenotypes', phenotype, '.', '.', n, p])

def meta():
	def first_append(d, i, name, count, percentage, sv, hap_score):
		if name not in d:
			d[name] = [sv, hap_score]
			for j in range(i):
				d[name].append('.')
				d[name].append('.')
		d[name].append(count)
		d[name].append(percentage)

	def second_append(d, i):
		for name in d:
			if len(d[name]) == 2 * i + 2:
				d[name].append('.')
				d[name].append('.')	

	args_.summary_list = [helper.check_file(x) for x in args_.summary_list]

	header1 = ['type', 'name', 'sv', 'hap_score']
	dicts = {}
	for i in range(len(args_.summary_list)):
		summary = args_.summary_list[i]
		prefix = os.path.basename(summary).replace('.result.summary.txt', '')
		header1.append(prefix + '_count')
		header1.append(prefix + '_percentage')
		with open(summary) as f:
			header2 = next(f).strip().split('\t')
			for line in f:
				fields = line.strip().split('\t')
				type = fields[header2.index('type')]
				name = fields[header2.index('name')]
				count = fields[header2.index('count')]
				percentage = fields[header2.index('percentage')]
				hap_score = fields[header2.index('hap_score')]
				sv = fields[header2.index('sv')]
				if type not in dicts:
					dicts[type] = {}
				first_append(dicts[type], i, name, count, percentage, sv, hap_score)
		for type in dicts:
			second_append(dicts[type], i)	
								
								
	with open(f'{args_.output_prefix}.stargazer-meta.txt', 'w') as f:
		f.write('\t'.join(header1) + '\n')
		
		def w(l):
			if len(l) != len(header1):
				raise ValueError()
			f.write('\t'.join([str(x) for x in l]) + '\n')		
		
		for type in dicts:
			if type == 'samples':
				for subtype in ['total', 'sv']:				
					w(['samples', subtype] + dicts[type][subtype])

			elif type == 'haps':
				for subtype in ['total', 'unique']:
					w(['haps', subtype] + dicts[type][subtype])
									
			elif type == 'alleles':				
				for allele in helper.sort_star_names(list(dicts[type])):
					w(['alleles', allele] + dicts[type][allele])
				
			elif type == 'phenotypes':
				for phenotype in sorted(list(dicts[type]), key = helper.PHENOTYPES.index):				
					w(['phenotypes', phenotype] + dicts[type][phenotype])
			
			else:
				for subtype in sorted(dicts[type]):
					w([type, subtype] + dicts[type][subtype])

def snp():
	if args_.vcf and args_.target_gene and args_.pair_list:
		pass
	else:
		raise ValueError('Missing one or more required arguments')
	f = open(args_.extended_prefix + '.txt', 'w')
	for pair in args_.pair_list:
		table = []
		finalized_vcf = helper.read_vcf_simple(args_.vcf)
		args_.samples = helper.vcf2samples(finalized_vcf)
		sample = args_.samples[pair.split('/')[0]]
		star = args_.star_dict[pair.split('/')[1]]		
		f.write("<sample={},star={}>\n".format(sample.name, star.name))
		header = ['pos', 'wt', 'var', 'type', 'so', 'impact', 'effect', 'hap1_allele', 'hap2_allele', 'gt', 'hap1_ad', 'hap2_ad', 'hap1_af', 'hap2_af']
		f.write('\t'.join(header) + '\n')

		def get_fields(snp, type):
			hap1_allele = snp.var if snp in sample.hap[0].obs else snp.wt
			hap2_allele = snp.var if snp in sample.hap[1].obs else snp.wt
			hap1_gt = '0' if hap1_allele == snp.wt else '1' if hap1_allele == snp.var else '2'
			hap2_gt = '0' if hap2_allele == snp.wt else '1' if hap2_allele == snp.var else '2'
			hap1_ad = str([x for x in sample.hap[0].obs if x.pos == snp.pos][0].ad) if snp.pos in [x.pos for x in sample.hap[0].obs] else '0'
			hap2_ad = str([x for x in sample.hap[1].obs if x.pos == snp.pos][0].ad) if snp.pos in [x.pos for x in sample.hap[1].obs] else '0'
			hap1_af = '{:.2f}'.format([x for x in sample.hap[0].obs if x.pos == snp.pos][0].af) if snp.pos in [x.pos for x in sample.hap[0].obs] else '0'
			hap2_af = '{:.2f}'.format([x for x in sample.hap[1].obs if x.pos == snp.pos][0].af) if snp.pos in [x.pos for x in sample.hap[1].obs] else '0'
			return [snp.pos, snp.wt, snp.var, type, snp.so, snp.impact, snp.effect, hap1_allele, hap2_allele, '{}|{}'.format(hap1_gt, hap2_gt), hap1_ad, hap2_ad, hap1_af, hap2_af]

		for snp in star.core:
			table.append(get_fields(snp, 'core'))
		for snp in star.tag:
			table.append(get_fields(snp, 'tag'))
		for fields in sorted(table, key = lambda x: int(x[0])):
			f.write('\t'.join(fields) + '\n')
		f.write('\n')
	f.close()

def plotcov():
	exit_code = subprocess.call(['Rscript', args_.program_dir + '/plotcov.R', args_.sdf, args_.output_prefix])

def changepoint():
	exit_code = subprocess.call(['Rscript', args_.program_dir + '/changepoint.R', args_.extended_prefix, args_.rdata])

SUBTOOLS = {'summary': summary, 'meta': meta, 'snp': snp, 'plotcov': plotcov, 'changepoint': changepoint}

DESCRIPTION = f'''tool description:
  perform secondary analyses of genotype calling

getting help:
  stargazer.py view -h

available subtools:
  summary
    create project-level summary using output data from the 'genotype' tool

  meta
    create a comparison chart by combining multiple reports from the 'report' subtool

  snp
    get variant data for given pairs of a sample and a star allele, by using 'finalized.vcf' generated by the 'genotype' tool

  changepoint
    detect changepoints (breakpoints) in copy number using 'sv.RData' generated by the 'genotype' tool

main usages:
  run the summary subtool
    stargazer.py view summary -o OUTPUT_PREFIX --genotype GENOTYPE

  run the meta subtool
    stargazer.py view meta -o OUTPUT_PREFIX --summary_list [SUMMARY [SUMMARY ...]]

  run the snp subtool
    stargazer.py view snp -o OUTPUT_PREFIX --pair_list [PAIR [PAIR ...]] --vcf VCF -t TARGET_GENE

  run the changepoint subtool
    stargazer.py view changepoint -o OUTPUT_PREFIX --rdata RDATA
'''

def run(args):
	global args_
	args_ = args

	SUBTOOLS[args_.subtool]()
					
