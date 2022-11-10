import helper

args_ = None

DESCRIPTION = f'''tool description:
  report genotype-based prescribing recommendations *testing only*

getting help:
  stargazer.py report -h

main usages:
  stargazer.py genotype -o OUTPUT_PREFIX --genotype GENOTYPE
'''

def get_meta():
	assessed_genes = len([x for x in args_.result_dict if args_.result_dict[x]['status'] != '-'])
	typed_genes = [args_.result_dict[x]['status'] for x in args_.result_dict].count('g')

	return '''<p>
  Sample ID: {id}<br />
  Date: {report_date}<br />
  Stargazer version: {version}<br />
  Genes examined: {assessed_genes}/{target_genes}<br />
  Genotypes called: {typed_genes}/{assessed_genes}<br />
</p>'''.format(id = args_.output_prefix, report_date = args_.date, version = args_.stargazer_version, assessed_genes = assessed_genes, target_genes = len(args_.target_genes), typed_genes = typed_genes)

def get_intro():
	return '''<p>
  Thank you for choosing Stargazer! Stargazer is a bioinformatiscs tool for predicting how a person's DNA affects their response to hundreds of medications. Stargazer does this by accurately calling star alleles (haplotypes) in pharmacogenetic (PGx) genes, which are defined by single-nucleotide variants (SNVs), small insertion-deletions (indels), and/or large structural variants (SVs). Once identified, these star alleles can be translated to an activity score (AS), which is in turn used to predict the person's drug response. Stargazer can utilize genomic data from various sources including next-generation sequencing (NGS) and single nucleotide polymorphism (SNP) array. For NGS data, Stargazer supports whole genome sequencing (WGS) as well as targeted sequencing such as whole exome sequencing (WES). For more details please visit the Stargazer website (https://stargazer.gs.washington.edu/stargazerweb/).<br />
  <br />
  This report includes hundreds of gene/drug pairs (e.g., <i>CYP2D6</i>/codeine) with accompanying levels of evidence for changing drug choice and dosing decisions. These pairs are described by the Clinical Pharmacogenetics Implementation Consortium (CPIC), which is an international assoication whose primary goal is to facilitate the use of PGx tests for patient care. Most importantly, CPIC provides detailed guidelines for helping clinicians understand how available genetic test results should be used to optimize drug therapy. As of {cpic_date}, there are {cpic_pairs} gene/drug pairs listed in the CPIC website (https://cpicpgx.org). Finally, the Food and Drug Administration (FDA) provides additional guidance by requiring applicable PGx test information be included in the drug labeling. These FDA-approved drug labels are included in this report as well.<br />
  <br />
  Disclaimer: This report is still very much in development. Please do not use it other than for code testing. Thank you.
</p>'''.format(cpic_date = args_.cpic_date, cpic_pairs = len(args_.cpic_list))

def get_overview_table():
	table = '''  <tr>
    <th style="width: 10%;">No.</th>
    <th style="width: 10%;">Gene</th>
    <th style="width: 20%;">Genotype</th>
    <th style="width: 10%;">Total AS</th>
    <th style="width: 30%;">Predicted Phenotype</th>
    <th style="width: 10%;">Drugs</th>
    <th style="width: 10%;">Guidelines</th>
  </tr>
'''

	for i, gene in enumerate(args_.target_genes, 1):
		color = 'black' if any([x in args_.result_dict[gene]['phenotype'] for x in ['normal', 'unknown', '-']]) else 'red'
		table += '''  <tr style="color: {color};">
    <td>{i}</td>
    <td><i>{gene}</i></td>
    <td><i>{genotype}</i></td>
    <td>{score}</td>
    <td>{phenotype}</td>
    <td>{drugs}</td>
    <td>{guidelines}</td>
  </tr>
'''.format(color = color, i = i, gene = gene.upper(), genotype = args_.result_dict[gene]['hap1_main'] + '/' + args_.result_dict[gene]['hap2_main'], score = args_.result_dict[gene]['dip_score'], phenotype = args_.result_dict[gene]['phenotype'], drugs = len([x['CPIC Level'] for x in args_.cpic_list if x['Gene'] == gene.upper()]), guidelines = len([x for x in args_.cpic_list if x['Gene'] == gene.upper() and x['Guideline'] != '-']))

	return '<table>\n' + table + '</table>'

def get_genotype_table():
	table = '''  <tr>
  	<th style="width: 5%;">No.</th>
  	<th style="width: 5%;">Gene</th>
    <th style="width: 20%;">Star Allele</th>
    <th style="width: 10%;">AS</th>
    <th style="width: 50%;">SNVs/Indels</th>
    <th style="width: 10%;">SVs</th>
  </tr>
'''

	for i, gene in enumerate(args_.target_genes, 1):
		hap1_main_core = '<br />'.join(args_.result_dict[gene]['hap1_main_core'].split(','))
		hap2_main_core = '<br />'.join(args_.result_dict[gene]['hap2_main_core'].split(','))
		table += '''  <tr>
    <td rowspan="2">{i}</td>
    <td rowspan="2"><i>{gene}</i></td>
    <td><i>{hap1_main}</i></td>
    <td>{hap1_score}</td>
    <td>{hap1_main_core}</td>
    <td>{hap1_sv}</td>
  </tr>
  <tr>
    <td><i>{hap2_main}</i></td>
    <td>{hap2_score}</td>
    <td>{hap2_main_core}</td>
    <td>{hap2_sv}</td>	
  </tr>
'''.format(i = i, gene = gene.upper(), hap1_main = args_.result_dict[gene]['hap1_main'], hap1_score = args_.result_dict[gene]['hap1_score'], hap1_main_core = hap1_main_core, hap1_sv = args_.result_dict[gene]['hap1_sv'], hap2_main = args_.result_dict[gene]['hap2_main'], hap2_score = args_.result_dict[gene]['hap2_score'], hap2_main_core = hap2_main_core, hap2_sv = args_.result_dict[gene]['hap2_sv'])

	return '<table>\n' + table + '</table>'

def get_drug_table():
	table = '''  <tr>
    <th style="width: 5%;">No.</th>
    <th style="width: 15%;">Drug</th>
    <th style="width: 5%;">Gene</th>
    <th style="width: 5%;">Level</th>
    <th style="width: 10%;">FDA</th>
    <th style="width: 60%;">Guideline</th>
  </tr>
'''

	for i, pair in enumerate(sorted(args_.cpic_list, key = lambda x: (x['Drug'].lower(), x['Gene'])), 1):
		bold = 'bold' if pair['Gene'].lower() in args_.result_dict else 'normal'
		color = 'black' if bold == 'normal' or any([x in args_.result_dict[pair['Gene'].lower()]['phenotype'] for x in ['normal', 'unknown', '-']]) else 'red'
		table += '''  <tr style="font-weight: {bold}; color: {color};">
    <td>{i}</td>
    <td>{drug}</td>
    <td><i>{gene}</i></td>
    <td>{level}</td>
    <td>{fda}</td>
    <td>{guideline}</td>
  </tr>
'''.format(bold = bold, color = color, i = i, drug = pair['Drug'], gene = pair['Gene'], level = pair['CPIC Level'], fda = pair['PGx on FDA Label'], guideline = pair['Guideline'])

	return '<table>\n' + table + '</table>'

def run(args):
	global args_
	args_ = args

	args_.result_dict = {}
	args_.cpic_list = []
	args_.cpic_date = ''

	with open(args_.genotype) as f:
		header = next(f).strip().split('\t')
		for line in f:
			fields = line.strip().split('\t')
			gene = fields[0]
			args_.result_dict[gene] = dict(zip(header, fields))
	
	for gene in args_.target_genes:
		if gene not in args_.result_dict:
			args_.result_dict[gene] = dict(zip(header, ['-' for x in header]))

	with open(args_.program_dir + '/cpicPairs.csv') as f:
		args_.cpic_date = next(f).strip().replace('"', '').replace('Date last updated: ', '')
		header = next(f).strip().split(',')
		for line in f:
			fields = line.replace(', ', '; ').strip().split(',')
			gene, drug = fields[0], fields[1]
			args_.cpic_list.append(dict(zip(header, [x if x else '-' for x in fields])))

	message = '''<!DOCTYPE html>
<html>
<head>
<title>Stargazer Report</title>
<style>
* {{
  font-family: Arial, Helvetica, sans-serif;
}}

table {{
  border-collapse: collapse;
  width: 100%;
  font-size: 80%;
}}

th, td {{
  border: 1px solid black;
  padding: 4px;
}}
</style>
</head>
<body>

<h1>Stargazer Report</h1>
{meta}
<h2>Introduction</h2>
{intro}

<h2>Sections</h2>
<ul>
  <li>Overview</li>
  <li>Genotypes</li>
  <li>Drugs</li>
</ul>

<p style="page-break-before: always;">

<h2>Overview</h2>
<p>PGx genes whose genotype leads to altered phenotype are shown in <span style="color: red;">red</span>.</p>
{overview_table}

<p style="page-break-before: always;">

<h2>Genotypes</h2>
{genotype_table}

<p style="page-break-before: always;">

<h2>Drugs</h2>
<p>Gene/drug pairs are shown in <span style="font-weight: bold;">bold</span> if genotype is available, and in <span style="font-weight: bold; color: red;">red</span> if altered phenotype is predicted.</p>
{drug_table}

</body>
</html>
'''.format(meta = get_meta(), intro = get_intro(), overview_table = get_overview_table(), genotype_table = get_genotype_table(), drug_table = get_drug_table())

	with open(args_.output_prefix + '.stargazer-report.html', 'w') as f:
		f.write(message)