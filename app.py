#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: piyus
"""

import numpy as np
from flask import Flask, request, jsonify, render_template, session
from flask_session import Session
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
import logging 
from flask_login import current_user
#from flask_login import 
import uuid 
from flask import (
    Flask,
    request,
    redirect,
    url_for,
    render_template,
    send_from_directory,
)
from werkzeug.utils import secure_filename


import matplotlib
matplotlib.use('agg')


def tbi(vcf):
    # pass
    os.chdir(UPLOAD_FOLDER)
    if vcf.endswith(".vcf.gz"):
    #print("zipped file ")
    
        print(vcf.split(".gz")[0])
        os.system("gunzip {vcf}".format(vcf=vcf))
        os.system("bgzip -c {vcf_un} > {vcf2}".format(vcf_un=vcf.split(".gz")[0] , vcf2 = vcf))
        os.system("tabix -p vcf {vcf}".format(vcf=vcf))

        

    elif vcf.endswith(".vcf"):
        print(vcf+".gz")
        print("Zipping using bgzip")
        os.system("bgzip -c {vcf} > {vcf2}".format(vcf=vcf,vcf2=vcf +".gz"))
        os.system("tabix -p vcf {vcf2}".format(vcf2 = vcf+".gz"))
        # print("index made")

     #vcf = str(vcf)+".gz"
     
    elif vcf.endswith(".txt.gz"):
        f = vcf.split(".txt.gz")[0]+".vcf"
        os.system("gunzip {vcf}".format(vcf=vcf))
        os.system("bcftools convert --tsv2vcf {v} -f ../hs37d5.fa -s {v}  -Ov -o {f}".format(v= vcf.split(".gz")[0], f = f))
        os.system("bgzip -c {f} > {vcf2}".format(f=f, vcf2=f +".gz"))
        os.system("tabix -p vcf {vcf2}".format(vcf2 = f+".gz"))
        

    else:
        #f = vcf.split("uploads/")[-1].split(".txt")[0]+".vcf"
        f = vcf.split(".txt")[0]+".vcf"
        print("F is ",f)
        os.system("bcftools convert --tsv2vcf {v} -f ../hs37d5.fa  -s {v}  -Ov -o {f}".format(v= vcf, f = f))
        print("ran!")
        os.system("bgzip -c {f} > {vcf2}".format(f=f, vcf2=f +".gz"))
        os.system("rm -r {f}".format(f=f))
        os.system("rm -r {v}".format(v=vcf))
        os.system("tabix -p vcf {vcf2}".format(vcf2 = f+".gz"))
        #os.system("mv {vcf2} ./uploads".format(vcf2=f +".gz"))
        #os.system("mv {tabix} ./uploads".format(tabix=f +".gz.tbi"))
  

    os.chdir(DIR_PATH)

dpe = []
final_features = []
sequence = ""
dp_features = []
dpe = []


UPLOAD_FOLDER = os.path.dirname(os.path.abspath(__file__)) + "/uploads/"


app = Flask(__name__, static_url_path="/static")
app.secret_key = "pimadi"
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_TYPE"] = "filesystem"
#sess = Session()
#sess.init_app(app)
#sess.init_app(app)

ALLOWED_EXTENSIONS = {"gz", "vcf"}

DIR_PATH = os.path.dirname(os.path.realpath(__file__))
app.config["MAX_CONTENT_LENGTH"] = 100 * 1024 * 1024


@app.route("/")
def home():
    #print(current_user.id)
    #print(session.sid)
    return render_template("home.html")


@app.route("/predictor/")
def predictor():
    return render_template("predictor.html")


@app.route("/about/")
def about():
    return render_template("contact.html")


@app.route("/predict", methods=["POST"])
def predict():
    # Add code to make .vcf.gz.tbi (V.IMP)
    try:
        ph_result = []
        recommend = None
        df = pypgx.load_recommendation_table()
        print("This will work")

        drug = request.form.get("drug")
        vcf = request.files.get("file")
        id = str(uuid.uuid4())
        session["username"] = id 
        print(request.files)
       # print(sess.sid)
        print(request.form)
        print(vcf)

        
        # saving the file
        filename = secure_filename(vcf.filename)

       # print(app.config["UPLOAD_FOLDER"])
        
        
        #path_upload = os.path.join("uploads/", filename)
        #print(path_upload)
        vcf.save(os.path.join(app.config["UPLOAD_FOLDER"], filename))
        print(os.getcwd())
        tbi("{name}".format(name=filename))
        # tbi("eee")
        name_ = filename.split(".")[0]+".vcf.gz"
        new_file_path =  "./uploads/{name}".format(name=name_)
        print("new_file path is:",new_file_path)
        cmd =  f"python3 run_pypgx.py -vcf {new_file_path} -j {id} -d {drug}"
        os.system(cmd)

        focus_df = df.loc[(df["Drug"] == str(drug).lower())]

        groups = focus_df.groupby(["Drug"])

        print(groups)

        genes = []
        for ind, row in focus_df.iterrows():
            genes.append(row["Gene1"])

        print(genes)

        unqiue_genes = list(set(genes))

        print(unqiue_genes)
        """
        if len(unqiue_genes) == 1:
            print("Running api")
            print(id)
            print(filename)
            print("going in",os.path.join(UPLOAD_FOLDER, name_))
            pypgx.api.pipeline.run_ngs_pipeline(
                str(unqiue_genes[0]).upper(),
                id,
                variants= new_file_path
            )
          

            print("ran!")
            print(id)
            """

        os.chdir(id)
            #print("unzip ")
        subprocess.run(
                shlex.split("unzip results.zip".format(dir=id))
            )  # change this and shutils
            #print("unzipped")

        file = [f for f in os.listdir("./") if f.startswith("tmp")]
            # file = list(filter(lambda x : x.startswith("tmp"), list(os.listdir("./"))))

        print(f"file : {file[0]}")

        df2 = pd.read_csv("{dir}/data.tsv".format(dir=file[0]), sep="\t")

        print(df2.head())

        genotype = str(df2["Genotype"])
        phenotype = df2["Phenotype"][0]
        print(phenotype)
        print(genotype)

        # try to avoid this
        os.chdir("../")

        # recommendation = str(df.loc[(df['Drug'] == str(drug).lower())])
        recommendation = focus_df[focus_df["Phenotype1"] == phenotype][
            "Recommendation"
        ]
        rec = recommendation.to_frame()

        print(rec.columns)

        df3 = rec.drop_duplicates(
            keep="first"
        )  # Removing duplicates and just keeping the first hit
        recommend = df3["Recommendation"].tolist()

        ph_result = []
        if phenotype.startswith("Normal Metabolizer") or phenotype.startswith(
            "Rapid Metabolizer"
        ):
            ph_result.append("Take it!")
        else:
            ph_result.append("Leave it!")

            # print(df3)

            # df3.to_csv('rec2.csv')

        # return df3.to_json()
    except Exception as e:
        logging.error(str(e))
        #return {"msg": str(e)}
        return render_template('error.html', message=str(e))
    
        

 
    
    
    os.system("rm -r {vcf}".format(vcf=new_file_path))
    os.system("rm -r {vcf_tbi}".format(vcf_tbi = new_file_path+".tbi"))
    os.system("rm -r {dir}".format(dir=id))
    

    return render_template("show.html", data=recommend, varchar=drug, diag=ph_result)


def allowed_file(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route("/submit", methods=["GET", "POST"])
def submit():
    if request.method == "POST":
        if "file" not in request.files:
            # print('No file attached in request')
            return render_template("predictor.html")
            return redirect(request.url)
        file = request.files["file"]
        if file.filename == "":
            # flash('No selected file')
            return render_template("predictor.html")
            return redirect(request.url)
        if file and allowed_file(file.filename):
           # session["username"] = file.filename
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config["UPLOAD_FOLDER"], filename))
            data, varchar = predict(
                os.path.join(app.config["UPLOAD_FOLDER"], filename), filename
            )
            # return redirect(url_for('uploaded_file', filename=filename))
            session.pop('username', None)
    # return render_template('predictor.html')
    # return render_template("show.html", data=data, varchar=varchar)


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
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


if __name__ == "__main__":
    app.run(host="0.0.0.0",debug=True)
