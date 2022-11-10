import helper
import statistics

def sdf2gdf():
	header = ['Locus', 'Total_Depth', 'Average_Depth_sample'] + ['Depth_for_' + x for x in args_.sample_list]

	# peek into the first line
	with open(args_.sdf) as f:
		fields = next(f).strip().split('\t')
		if len(fields) != len(header) - 1:
			raise ValueError('Input SDF and sample list different in length')

	with open(args_.output_prefix + '.gdf', 'w') as f1:
		f1.write('\t'.join(header) + '\n')
		with open(args_.sdf) as f2:
			for line in f2:
				fields = line.strip().split('\t')
				chr = fields[0].replace('chr', '')
				pos = fields[1]
				locus = chr + ':' + pos
				data = [int(x) for x in fields[2:]]
				total = sum(data)
				avg = round(statistics.mean(data), 2)
				new_fields = [locus, total, avg] + data
				f1.write('\t'.join([str(x) for x in new_fields]) + '\n')

def define():
	with open(args_.table) as f:
		s = 14 # first star allele index
		header = next(f).strip().split('\t')
		dat = [[] for x in header[s:]]
		for line in f:
			fields = line.strip().split('\t')
			pos = fields[5]
			var = fields[8]
			wt = fields[9]
			for i in range(s, len(fields)):
				if fields[i] == '1':
					dat[i - s].append('{}:{}>{}'.format(pos, wt, var))

	with open(args_.output_prefix + '.stargazer-view.txt', 'w') as f:
		for i in range(len(dat)):
			f.write(header[s + i] + '\t' + ','.join(dat[i]) + '\n')

def slice():
	sliced_vcf = helper.read_vcf_region(args_.vcf, args_.region)
	helper.write_vcf(sliced_vcf, args_.output_prefix + '.vcf')

SUBTOOLS = {'sdf2gdf': sdf2gdf, 'define': define, 'slice': slice}

DESCRIPTION = f'''tool description:
  create various files necessary for running Stargazer

getting help:
  stargazer.py setup -h

available subtools:
  sdf2gdf
    create gdf file from sdf file
    
  define
    define star alleles using variants from the SNP table
    
  slice
    create sliced VCF file

main usages:
  run the sdf2gdf subtool
    stargazer.py setup sdf2gdf -o OUTPUT_PREFIX --sample_list [SAMPLE [SAMPLE ...]] --sdf SDF
  
  run the define subtool
    stargazer.py setup define -o OUTPUT_PREFIX --table TABLE
  
  run the slice subtool
    stargazer.py setup slice -o OUTPUT_PREFIX --region REGION --VCF VCF
'''

def run(args):
	global args_
	args_ = args

	SUBTOOLS[args_.subtool]()