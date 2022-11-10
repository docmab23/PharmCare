import argparse, datetime, os, sys, timeit, types
import helper, genotype, setup, view, pipeline, report

ABBREVIATIONS = 'abbreviations and acronyms:\n  1KGP, The 1000 Genomes Project; AF, allele fraction; AS, activity score; CN, copy number; CNSR, copy number-stable region; DOI, digital object identifier; DPSV, The Database of Pharmacogenomic Structural Variants; GATK, The Genome Analysis Toolkit; GDF, GATK-DepthOfCoverage format; indels, insertion-deletion variants; NGS, next-generation sequencing; PGx, pharmacogenomics; SDF, SAMtools depth format; SGE, Sun Grid Engine; SNP, single nucleotide polymorphism; SNV, single nucleotide variant; SV, structural variant; TS, targeted sequencing; VCF, variant call format; WGS, whole genome sequencing'
TOOLS = {'genotype': genotype, 'setup': setup, 'view': view, 'pipeline': pipeline, 'report': report}
DESCRIPTION = f'''available tools:
  genotype
    call star alleles in target gene from genomic data

  setup
    create various files necessary for running Stargazer

  view
    perform secondary analyses of genotype calling

  pipeline
    provide end-to-end solutions for genotyping pipeline

  report
    report genotype-based prescribing recommendations *testing only*

getting help:
  stargazer.py -h
  stargazer.py tool -h
'''

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\n***********************************************************************\n%s: error: %s\n***********************************************************************\n' % (self.prog, message))

def main():
	# get the start timestamp
	start_time = timeit.default_timer()

	# create the top-level parser
	parser = ArgumentParser(formatter_class = argparse.RawTextHelpFormatter, description = DESCRIPTION, epilog = f'{ABBREVIATIONS}', usage = argparse.SUPPRESS)
	subparsers = parser.add_subparsers(dest = 'tool', metavar = 'tool', help = 'tool name')
	subparsers.required = True

	# create the parent parser
	parent_parser = argparse.ArgumentParser(add_help = False)
	parent_parser.add_argument('--output_dir', type = helper.check_file, help = 'output files directory (default: current working directory)')
	parent_parser.add_argument('-o', '--output_prefix', required = True, help = 'output filename prefix')

	# create the parser for the "genotype" tool
	genotype_parser = subparsers.add_parser('genotype', parents = [parent_parser], formatter_class = argparse.RawTextHelpFormatter, epilog = f'{ABBREVIATIONS}', description = genotype.DESCRIPTION, usage = argparse.SUPPRESS)
	genotype_parser.add_argument('-c', '--control_gene', type = str.lower, choices = helper.get_control_genes(), metavar = 'CONTROL_GENE', help = 'control gene name')
	genotype_parser.add_argument('--control_type', type = str.lower, choices = ['known', 'cnsr', 'custom'], default = 'known', help = "type of the control locus for intra-sample normalization of read depth (choose from 'known', 'cnsr', 'custom') (default: known)", metavar = 'CONTROL_TYPE')
	genotype_parser.add_argument('-d', '--data_type', required = True, type = str.lower, choices = ['ts', 'wgs', 'chip'], metavar = 'DATA_TYPE', help = "type of input data (choose from 'wgs', 'ts', 'chip')")
	genotype_parser.add_argument('--detail', action = 'store_true', help = 'output the detailed version of the CN and AF plots')
	genotype_parser.add_argument('--gdf', type = helper.check_file, help = 'input GDF file')
	genotype_parser.add_argument('--impute', action = 'store_true', help = 'allow imputation of ungenotyped markers during statistical phasing')	
	genotype_parser.add_argument('--ref_vcf', help = 'reference haplotype panel file (.vcf or .gz) (default: 1KGP phase 3)')
	genotype_parser.add_argument('--region', help = 'genomic region (format: chr:start-end)')	
	genotype_parser.add_argument('--sample_list', nargs = '*', metavar = 'SAMPLE', help = 'space-delimited list of sample IDs')	
	genotype_parser.add_argument('-t', '--target_gene', required = True, choices = helper.get_target_genes(), metavar = 'TARGET_GENE', type = str.lower, help = 'target gene name')
	genotype_parser.add_argument('--vcf', required = True, type = helper.check_file, help = 'input VCF file (.vcf or .vcf.gz)')

	# create the parser for the "setup" tool
	setup_parser = subparsers.add_parser('setup', parents = [parent_parser], formatter_class = argparse.RawTextHelpFormatter, usage = argparse.SUPPRESS, description = setup.DESCRIPTION)
	setup_parser.add_argument('subtool', type = str.lower, choices = setup.SUBTOOLS, metavar = 'subtool', help = 'subtool name')
	setup_parser.add_argument('--region', help = 'genomic region (format: chr:start-end)')
	setup_parser.add_argument('--sample_list', nargs = '*', metavar = 'SAMPLE', help = 'space-delimited list of sample IDs')
	setup_parser.add_argument('--sdf', type = helper.check_file, help = 'input SDF file')	
	setup_parser.add_argument('--table', type = helper.check_file, help = "modified SNP table ('snp_table.txt') with additional columns, each representing a star allele to be defined (put '1' if the variant is used, '0' if otherwise)")
	setup_parser.add_argument('--vcf', type = helper.check_file, help = 'input VCF file (.vcf or .vcf.gz)')

	# create the parser for the "view" tool
	view_parser = subparsers.add_parser('view', parents = [parent_parser], formatter_class = argparse.RawTextHelpFormatter, usage = argparse.SUPPRESS, description = view.DESCRIPTION)
	view_parser.add_argument('subtool', type = str.lower, choices = view.SUBTOOLS, metavar = 'subtool', help = 'subtool name')
	view_parser.add_argument('--genotype', type = helper.check_file, help = 'genotype call file')
	view_parser.add_argument('--pair_list', nargs = '*', metavar = 'PAIR', help = 'space-delimited list of pairs of a sample and a star allele (format: sample/star allele)')
	view_parser.add_argument('--rdata', type = helper.check_file, help = 'input RData file')	
	view_parser.add_argument('--sdf', type = helper.check_file, help = 'input SDF file')
	view_parser.add_argument('--summary_list', nargs = '*', metavar = 'SUMMARY', help = 'space-delimited list of summary files')
	view_parser.add_argument('-t', '--target_gene', choices = helper.get_target_genes(), metavar = 'TARGET_GENE', type = str.lower, help = 'target gene name')	
	view_parser.add_argument('--vcf', type = helper.check_file, help = 'input VCF file (.vcf or .vcf.gz)')
		
	# create the parser for the "pipeline" tool
	pipeline_parser = subparsers.add_parser('pipeline', parents = [parent_parser], formatter_class = argparse.RawTextHelpFormatter, usage = argparse.SUPPRESS, epilog = f'{ABBREVIATIONS}', description = pipeline.DESCRIPTION)
	pipeline_parser.add_argument('subtool', type = str.lower, choices = pipeline.SUBTOOLS, metavar = 'subtool', help = 'subtool name')
	pipeline_parser.add_argument('--assembly', required = True, type = helper.check_file, help = 'reference genome assembly file (.fa or .fasta)')
	pipeline_parser.add_argument('--bam', type = helper.check_file, help = 'input BAM file')
	pipeline_parser.add_argument('--bam_list', nargs = '*', metavar = 'BAM', help = 'space-delimited list of BAM files')
	pipeline_parser.add_argument('--batch', type = int, default = 1, help = 'number of batches used for parallelization (default: 1)')	
	pipeline_parser.add_argument('-c', '--control_gene', required = True, help = 'control gene name')
	pipeline_parser.add_argument('-d', '--data_type', type = str.lower, choices = ['ts', 'wgs', 'chip'], metavar = 'DATA_TYPE', help = "type of input data (choose from 'wgs', 'ts', 'chip')")
	pipeline_parser.add_argument('--dbsnp', required = True, type = helper.check_file, help = 'dbSNP file (.vcf or .vcf.gz)')
	pipeline_parser.add_argument('--gatk', type = helper.check_file, help = 'GATK3 file (GenomeAnalysisTK.jar)')	
	pipeline_parser.add_argument('--gene_list', nargs = '*', default = [], metavar = 'GENE', help = 'space-delimited list of target genes (default: all)')	
	pipeline_parser.add_argument('--manifest', type = helper.check_file, help = 'manifest file')
	pipeline_parser.add_argument('--mapping', type = int, default = 1, help = 'minimum read mapping quality score used when computing read depth')
	pipeline_parser.add_argument('--sample_list', nargs = '*', metavar = 'SAMPLE', help = 'space-delimited list of sample IDs')
	pipeline_parser.add_argument('--sampling', action = 'store_true', help = 'subset the control locus when computing read depth to speed up processing')
	pipeline_parser.add_argument('-t', '--target_gene', choices = helper.get_target_genes(), metavar = 'TARGET_GENE', type = str.lower, help = 'target gene name')

	# create the parser for the "report" tool
	report_parser = subparsers.add_parser('report', parents = [parent_parser], formatter_class = argparse.RawTextHelpFormatter, usage = argparse.SUPPRESS, description = report.DESCRIPTION)
	report_parser.add_argument('--genotype', type = helper.check_file, help = 'genotype call file')

	# parse the arguments
	args = parser.parse_args()
	
	# update the namespace for variables shared by more than one tool
	args.program_dir = os.path.dirname(os.path.realpath(__file__))
	args.stargazer_version = args.program_dir.split('/')[-1]
	args.tool_info = f'stargazer-{args.tool}-{args.subtool}' if hasattr(args, 'subtool') else f'stargazer-{args.tool}'
	args.extended_prefix = f'{args.output_prefix}.{args.tool_info}'
	args.log = f'{args.extended_prefix}.log'
	args.date = datetime.datetime.now().strftime('%Y%m%d')
	args.time = datetime.datetime.now().strftime('%H:%M:%S')
	args.line_break = '+{}\n'.format('-' * 70)
	args.target_region = '' if not hasattr(args, 'target_gene') or not args.target_gene else helper.get_region(args.target_gene)
	args.control_region = helper.get_region(args.control_gene) if hasattr(args, 'control_gene') and args.control_gene else helper.get_cnsr(args.target_gene) if hasattr(args, 'control_type') and args.control_type == 'cnsr' else args.region.replace('chr', '') if hasattr(args, 'control_type') and args.control_type == 'custom' else ''	
	args.snp_list = helper.get_snp_list(args.target_gene) if hasattr(args, 'target_gene') and args.target_gene else '.'
	args.star_dict = helper.get_star_dict(args.target_gene, args.snp_list) if hasattr(args, 'target_gene') and args.target_gene else '.'
	args.samples = {}
	args.target_genes = helper.get_target_genes()
	args.gene_dict = helper.get_gene_dict()

	# move to the output directory
	if args.output_dir:
		os.chdir(args.output_dir)

	# write intro
	intro = '\n'.join([args.stargazer_version, 'Author: Seung-been "Steven" Lee', 'Enter "python3 stargazer.py --help" to view command line arguments', 'For more details, please visit https://stargazer.gs.washington.edu/stargazerweb/index.html', '', 'Stargazer tool used: ' + args.tool, 'Date: ' + args.date, 'Start time: ' + args.time, '', 'Command line: python3 ' + ' '.join(sys.argv), ''])
	print('\n' + intro)
	with open(args.log, 'w') as f:
		f.write(intro)

	# run the selected tool
	TOOLS[args.tool].run(args)

	# write outro
	stop_time = timeit.default_timer()
	elapsed_time = str(datetime.timedelta(seconds = (stop_time - start_time))).split('.')[0]
	outro = '\n'.join([f'Elapsed time: {elapsed_time}', 'Stargazer finished'])
	print(outro + '\n')
	with open(args.log, 'a') as f:
		f.write('\n')
		f.write(outro)

if __name__ == '__main__':
	main()
