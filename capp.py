#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: piyus
"""

import numpy as np
from flask import Flask, request, jsonify, render_template
import pickle
import pypgx 
from pypgx.api import utils
import pandas as pd
import subprocess
import argparse
import os 
import shlex
from werkzeug.utils import secure_filename
import os
from flask import Flask, request, redirect, url_for, render_template, send_from_directory
from werkzeug.utils import secure_filename


dpe=[]
final_features=[]
sequence=""
dp_features=[]
dpe=[]

UPLOAD_FOLDER = os.path.dirname(os.path.abspath(__file__)) + '/uploads/'


app = Flask(__name__, static_url_path="/static")
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

ALLOWED_EXTENSIONS = {'fasta', 'vcf'}

DIR_PATH = os.path.dirname(os.path.realpath(__file__))
app.config['MAX_CONTENT_LENGTH'] = 5 * 1024 * 1024

@app.route('/')
def home():
	return render_template('home.html')

@app.route('/predictor/')
def predictor():
    return render_template('predictor.html')

@app.route('/about/')
def about():
    return render_template('contact.html')





@app.route('/predict',methods=['POST'])
def predict():
	# Add code to make .vcf.gz.tbi (V.IMP)
	try:
		df = pypgx.load_recommendation_table()
		print("This will work")
		print(request.files)

		drug = request.form.get("drug")
		vcf = request.files.get('file')
		id = "dekhe"
		print(drug)

		# saving the file
		filename = secure_filename(vcf.filename)
		print(app.config['UPLOAD_FOLDER'])
		vcf.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

		focus_df = df.loc[(df['Drug'] == str(drug).lower())]

		groups = (focus_df.groupby(["Drug"]))

		print(groups)

		genes=[]
		for ind , row in focus_df.iterrows():
			genes.append(row["Gene1"])

		print(genes)

		unqiue_genes = list(set(genes))

		print(unqiue_genes)

		if len(unqiue_genes)==1:
			print("Running api")
			print(id)
			print(filename)
			pypgx.api.pipeline.run_ngs_pipeline(str(unqiue_genes[0]).upper() ,id , variants=os.path.join(UPLOAD_FOLDER, filename))

			print(id)

			os.chdir(id)
			print("unzip ")
			subprocess.run(shlex.split("unzip results.zip".format(dir=id))) #change this and shutils
			print("unzipped")
			
			file = [f for f in os.listdir("./")  if f.startswith("tmp") ]
			#file = list(filter(lambda x : x.startswith("tmp"), list(os.listdir("./"))))

			print(f"file : {file[0]}")
			
			df2 = pd.read_csv("{dir}/data.tsv".format(dir=file[0]),sep='\t')
			
			print(df2.head())

			genotype = str(df2["Genotype"])
			phenotype = (df2["Phenotype"][0])
			print(phenotype)
			print(genotype)

			# try to avoid this			
			os.chdir("../")

			

			#recommendation = str(df.loc[(df['Drug'] == str(drug).lower())])
			recommendation = focus_df[focus_df["Phenotype1"]==phenotype]["Recommendation"]
			rec=recommendation.to_frame()
			
			print(rec.columns)
			
			df3 = rec.drop_duplicates(keep='first') #Removing duplicates and just keeping the first hit
			recommend = df3["Recommendation"].tolist()

			ph_result=[]
			if phenotype.startswith("Normal Metabolizer") or phenotype.startswith("Rapid Metabolizer"):
				ph_result.append("Take it!")
			else:
				ph_result.append("Leave it!")



       		
     		
       		   

			#print(df3)
			
			
			#df3.to_csv('rec2.csv')
		
		#return df3.to_json()
	

	except Exception as e:
		print(str(e))
		return {"msg" : str(e)}
	
	return render_template("show.html", data=recommend, varchar=drug, diag=ph_result)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/submit', methods=['GET', 'POST'])
def submit():
    if request.method == 'POST':
        if 'file' not in request.files:
            #print('No file attached in request')
            return render_template('predictor.html')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            #flash('No selected file')
            return render_template('predictor.html')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            data,varchar=predict(os.path.join(app.config['UPLOAD_FOLDER'], filename), filename)
            #return redirect(url_for('uploaded_file', filename=filename))
    #return render_template('predictor.html')
    #return render_template("show.html", data=data, varchar=varchar)

"""
def process_file(path, filename):
	drug =  request.form.get("drug")
	#vcf = args.vcf
	#id = args.jobid

	focus_df = df.loc[(df['Drug'] == str(drug).lower())]

	groups = (focus_df.groupby(["Drug"]))

	genes=[]
	for ind , row in focus_df.iterrows():
    	genes.append(row["Gene1"])


	unqiue_genes = list(set(genes))


	if len(unqiue_genes)==1:
    	pypgx.api.pipeline.run_ngs_pipeline(str(unqiue_genes[0]).upper() ,id , variants=vcf)
    #os.chdir("./{dir}".format(dir=id))
    #print("pipeline done !")
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

     
     return(str,kratikal)
     #return render_template("show.html", data=str, varchar=kratikal)

"""



def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS



if __name__ == "__main__":
    app.run(debug=True)

