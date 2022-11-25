import pypgx 
from pypgx.api import utils
import pandas as pd
import subprocess
import argparse
import os 
import shlex

df = pypgx.load_recommendation_table()

#df.to_csv("utils/recommend.csv")
parser = argparse.ArgumentParser(prog = 'pypgx')
#The ArgumentParser.add_argument() method attaches individual argument specifications to the parser. It supports positional arguments, options that accept values, and on/off flags:

#parser.add_argument('filename')           # positional argument
parser.add_argument('-vcf', '--vcf' , required=True)      # option that takes a value
parser.add_argument('-j', '--jobid', required=True)
parser.add_argument('-d', '--drug' , required=True)    

#genes =['cacna1s', 'cftr', 'cyp1a1', 'cyp1a2', 'cyp1b1', 'cyp2a6', 'cyp2a13', 'cyp2b6', 'cyp2c8', 'cyp2c9', 'cyp2c19', 'cyp2d6', 'cyp2e1', 'cyp2f1', 'cyp2j2', 'cyp2r1', 'cyp2s1', 'cyp2w1', 'cyp3a4', 'cyp3a5', 'cyp3a7', 'cyp3a43']
args = parser.parse_args()

drug = args.drug 
vcf = args.vcf
id = args.jobid

focus_df = df.loc[(df['Drug'] == str(drug).lower())]

groups = (focus_df.groupby(["Drug"]))

genes=[]
for ind , row in focus_df.iterrows():
    genes.append(row["Gene1"])


unqiue_genes = list(set(genes))
print(unqiue_genes)


if len(unqiue_genes)==1:
    pypgx.api.pipeline.run_ngs_pipeline(str(unqiue_genes[0]).upper() ,id , variants=vcf)
    #os.chdir("./{dir}".format(dir=id))
    print("pipeline done !")
    os.chdir(id)
    subprocess.run(shlex.split("unzip results.zip".format(dir=id)))
    print("unzipped")


    file = [f for f in os.listdir("./")  if f.startswith("tmp") ]

    print(file[0])
    df2 = pd.read_csv("{dir}/data.tsv".format(dir=file[0]),sep='\t')
    
    genotype = str(df2["Genotype"])
    phenotype = (df2["Phenotype"][0])
    print(phenotype)

    os.chdir("../")
    

    #recommendation = str(df.loc[(df['Drug'] == str(drug).lower())])
    recommendation = focus_df[focus_df["Phenotype1"]==phenotype]["Recommendation"]
    rec=recommendation.to_frame()
    print(rec.columns)
    df3 = rec.drop_duplicates(keep='first') #Removing duplicates and just keeping the first hit
    print(df3)
    df3.to_csv('rec2.csv')












    


#print(groups)

#print(pypgx.predict_phenotype('CFTR', '*1', '*5'))

#print(pypgx.get_recommendation('sertraline', 'CYP2C19', 'Ultrarapid Metabolizer'))
#


#pypgx.api.pipeline.run_ngs_pipeline("CYP1A1" , "new_fi", variants="test.vcf.gz")
#df = pypgx.load_cpic_table()
#df.to_csv("utils/cpic_table.csv")
