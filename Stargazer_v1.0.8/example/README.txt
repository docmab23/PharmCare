# Get the path to the Stargazer directory

x=/Users/seungbeenlee/Desktop/Stargazer_v1.0.8

# Move to the example directory

cd $x/example

# Display all of the Stargazer tools

python3 $x/stargazer.py --help

# Display tool-specific arguments for, say, the 'genotype' tool

python3 $x/stargazer.py genotype --help

##############################################################################
# Genotype analysis for WGS data (GeT-RM; N=70) from Lee et al. (2019)       #
##############################################################################

# Below uses VDR as the control gene for copy number (CN) analysis

python3 $x/stargazer.py genotype \
  -o getrm-cyp2d6-vdr \
  -d wgs \
  -t cyp2d6 \
  -c vdr \
  --vcf getrm-cyp2d6-vdr.joint.filtered.vcf \
  --gdf getrm-cyp2d6-vdr.gdf

# It's recommended to check CN results using multiple control genes
# Additionally, you may use the copy number-stable region (CNSR) as control

python3 $x/stargazer.py genotype \
  -o getrm-cyp2d6-cnsr \
  -d wgs \
  -t cyp2d6 \
  --control_type cnsr \
  --vcf getrm-cyp2d6-vdr.joint.filtered.vcf \
  --gdf getrm-cyp2d6-vdr.gdf

# Finally, you can provide a custom region as control

python3 $x/stargazer.py genotype \
  -o getrm-cyp2d6-custom \
  -d wgs \
  -t cyp2d6 \
  --control_type custom \
  --region chr22:42546883-42551883 \
  --vcf getrm-cyp2d6-vdr.joint.filtered.vcf \
  --gdf getrm-cyp2d6-vdr.gdf

##############################################################################
# Genotype analysis for TS data (PGRNseq; N=96) from Lee et al. (2018)       #
##############################################################################

# Unlike WGS data, TS data requires inter-sample normalization for CN analysis
# Below uses the population mean during inter-sample normalization

python3 $x/stargazer.py genotype \
  -o hapmap-cyp2d6-vdr \
  -d ts \
  -t cyp2d6 \
  -c vdr \
  --vcf hapmap-cyp2d6-vdr.joint.filtered.vcf \
  --gdf hapmap-cyp2d6-vdr.gdf
  
# For CN analysis, you may indicate known reference samples without SV
# Below will use the mean of indicated samples instead of the population mean

python3 $x/stargazer.py genotype \
  -o hapmap-cyp2d6-vdr-ref \
  -d ts \
  -t cyp2d6 \
  -c vdr \
  --sample_list 133419 133420 133421 133423 133425 \
  --vcf hapmap-cyp2d6-vdr.joint.filtered.vcf \
  --gdf hapmap-cyp2d6-vdr.gdf

##############################################################################
# Genotype analysis with a VCF file only (WGS data)                          #
##############################################################################

python3 $x/stargazer.py genotype \
  -o getrm-cyp3a5-vcfonly \
  -d wgs \
  -t cyp3a5 \
  --vcf getrm-cyp3a5-vdr.joint.filtered.vcf

##############################################################################
# Genotype analysis with a VCF file only (SNP array data)                    #
##############################################################################

python3 /Users/seungbeenlee/Desktop/Stargazer_v1.0.8/stargazer.py genotype \
  -o rok-cyp3a5-chip \
  -d chip \
  -t cyp3a5 \
  --vcf rok-cyp3a5.vcf
