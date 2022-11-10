import os 
import sys 
import argparse
import subprocess
import shlex 

stargazer_dir = "Stargazer_v1.0.8"
res_dir = "stargazer_results"
parser = argparse.ArgumentParser(
                    prog = 'stargazer')
#The ArgumentParser.add_argument() method attaches individual argument specifications to the parser. It supports positional arguments, options that accept values, and on/off flags:

#parser.add_argument('filename')           # positional argument
parser.add_argument('-vcf', '--vcf' , required=True)      # option that takes a value
#parser.add_argument('-g', '--gene', required=True)
genes =['cacna1s', 'cftr', 'cyp1a1', 'cyp1a2', 'cyp1b1', 'cyp2a6', 'cyp2a13', 'cyp2b6', 'cyp2c8', 'cyp2c9', 'cyp2c19', 'cyp2d6', 'cyp2e1', 'cyp2f1', 'cyp2j2', 'cyp2r1', 'cyp2s1', 'cyp2w1', 'cyp3a4', 'cyp3a5', 'cyp3a7', 'cyp3a43']
args = parser.parse_args()

stargazer_file = stargazer_dir + "/stargazer.py"
print(stargazer_file)

#os.chdir(stargazer_dir)

#print("changed")
#print(str(stargazer_file))
print(args.vcf )

if not os.path.exists(res_dir):
    os.mkdir(res_dir)

#print(shlex.split("python stargazer.py genotype -o {r} -d wgs -t {t} --vcf {v}".format(s=stargazer_file, r=res_dir, t=args.gene, v=args.vcf)))


print(os.getcwd())

for ind , gene in enumerate(genes):
    print(gene)
    subprocess.run(shlex.split("python {s} genotype --output_dir {r} -o {ind} -d wgs -t {t} --vcf {v}".format(s=stargazer_file, r=res_dir, t=gene, v=args.vcf,ind=ind)))