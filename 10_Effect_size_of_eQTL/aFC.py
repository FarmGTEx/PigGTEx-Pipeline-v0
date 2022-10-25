#!/usr/bin/env python3

import argparse
import copy
import gzip
import math
import numpy as np
import os
import pandas as pd
import pysam
import scikits.bootstrap as boot
import statsmodels.api as sm
import statsmodels.formula.api as smf
import subprocess
import tempfile
import time
import warnings


def main():
	parser = argparse.ArgumentParser()
	# REQUIRED
	parser.add_argument("--vcf", required=True, help="Genotype VCF")
	parser.add_argument("--pheno", required=True, help="Phenotype file")
	parser.add_argument("--qtl", required=True, help="File containing QTL to calculate allelic fold change for. Should contain tab separated columns 'pid' with phenotype (gene) IDs and 'sid' with SNP IDs. Optionally can include the columns 'sid_chr' and 'sid_pos', which will facilitate tabix retrieval of genotypes, greatly reducing runtime.")
	parser.add_argument("--geno", required=False, default="GT", help="Which field in VCF to use as the genotype. By default 'GT' = genotype. Setting to 'DS' will use dosage rounded to the nearest integer (IE 1.75 = 2 = 1|1).")
	parser.add_argument("--chr", type=str, help="Limit to a specific chromosome.")
	parser.add_argument("--log_xform", type=int, required=True, help="The data has been log transformed (1/0). If so, please set --log_base.")
	parser.add_argument("--output", "--o", required=True, help="Output file")

	# OPTIONAL
	parser.add_argument("--cov", help="Covariates file")
	parser.add_argument("--matrix_o", help="Output the raw data matrix used to calculate aFC for each QTL into the specific folder.")
	parser.add_argument("--boot", default=100, type=int, help="Number of bootstraps to perform for effect size confidence interval. Can be set to 0 to skip confidence interval calculation, which will greatly reduce runtimes.")
	parser.add_argument("--ecap", default=math.log(100,2), type=float, help="Absolute aFC cap in log2.")
	parser.add_argument("--log_base", default=2, type=int, help="Base of log applied to data. If other than 2, data will be converted to log2.")
	parser.add_argument("--min_samps", default=2, type=int, help="Minimum number of samples with genotype data required to calculate effect size, default = 2.")
	parser.add_argument("--min_alleles", default=1, type=int, help="Minimum observations of each allele in data to calculate aFC, default = 1.")
	parser.add_argument("--count_o", default=0, type=int, help="Output the observed allele counts, default = 0.")

	# disable warnings
	warnings.filterwarnings("ignore")

	global args

	args = parser.parse_args()

	version = "0.3"
	print("")
	print("########################################################")
	print("		 Welcome to aFC v%s"%(version))
	print("  Authors: Pejman Mohammadi (pmohammadi@nygenome.org),\n	   Stephane Castel (scastel@nygenome.org)")
	print("########################################################")
	print("")

	print("RUN SETTINGS")
	print("     Genotype VCF: %s"%(args.vcf))
	print("     Phenotype File: %s"%(args.pheno))
	if args.cov != None:
		print("     Covariate File: %s"%(args.cov))
	print("     QTL File: %s"%(args.qtl))
	print("     Genotype Field: %s"%(args.geno))
	print("     Log Transformed: %d"%(args.log_xform))
	if args.log_xform == 1:
		print("     Log Base: %d"%(args.log_base))
	print("     Boostraps: %d"%(args.boot))
	print("     Minimum number of samples with genotype data: %d"%(args.min_samps))
	print("     Minimum number of allele observations: %d"%(args.min_alleles))
	if args.chr != None:
		print("     Chromosome: %s"%(args.chr))

	print("")

	if args.log_xform == 1:
		print("!! PLEASE ENSURE THAT YOUR DATA HAVE BEEN LOG TRANSFORMED WITH A BASE OF %d !!"%args.log_base)
	elif args.log_xform == 0:
		print("!! PLEASE ENSURE THAT YOUR DATA HAVE NOT BEEN LOG TRANSFORMED !!")

	print("")

	start_time = time.time()

	# get sample - column map from VCF
	print("1. Loading VCF...")
	vcf_map = sample_column_map(args.vcf)
	tabix_vcf = pysam.Tabixfile(args.vcf,"r")

	global df_cov
	if args.cov != None:
		print("1b. Loading covariates...")
		df_cov = pd.read_csv(args.cov, sep="\t", index_col=False)
		if "ID" in df_cov.columns:
			cov_id_col = "ID"
		elif "id" in df_cov.columns:
			cov_id_col = "id"
		else:
			print("First column was used.")
			# print(df_cov.columns)
			df_cov.rename(columns={'Unnamed: 0':'ID'},inplace=True)
			# print(df_cov.columns)
			cov_id_col = "ID"
	else:
		df_cov = pd.DataFrame(columns=['ID'])
		cov_id_col = "ID"

	#2 get sample - column map from phenotype file
	print("2. Loading phenotype data...")
	pheno_map = sample_column_map(args.pheno, line_key="#", start_col=4)
	tabix_pheno = pysam.Tabixfile(args.pheno, "r")

	# 3 load fastQTL results
	print("3. Loading fastQTL results...")
	df_qtl = pd.read_csv(args.qtl, sep="\t", index_col=False)

	print("4. Retrieving eSNP positions...")
	set_esnp = set(df_qtl['variant_id'].tolist())
	dict_esnp = {}

	if "sid_chr" in df_qtl.columns and "sid_pos" in df_qtl.columns:
		# eSNP positions are specified in file
		for index, row in df_qtl.iterrows():
			if args.chr == None or str(row['sid_chr']) == args.chr:
				dict_esnp[row['variant_id']] = [row['sid_chr'],int(row['sid_pos'])]
	else:
		# retrieve SNP positions from the VCF (since these are not included in the fastQTL output)
		print("     unpacking VCF...")
		# retrieve the SNP positions from the VCF
		tfile = tempfile.NamedTemporaryFile(delete=False)
		vcf_in = tfile.name
		tfile.close()
		if args.chr != None:
			# retrieve only genotypes from the desired chromosome
			error = subprocess.call("tabix "+args.vcf+" "+args.chr+": | cut -f 1-3 > "+vcf_in, shell=True)
			if error != 0:
				print("     ERROR loading retrieving genotype data. Ensure tabix index exists and is current.")
				quit()
		else:
			error = subprocess.call("gunzip -c "+args.vcf+" | cut -f 1-3 > "+vcf_in, shell=True)
			if error != 0:
				print("     ERROR loading retrieving genotype data.")
				quit()

		stream_in = open(vcf_in, "r")

		current_chr = ""

		for line in stream_in:
			if line[0:1] != "#":
				#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
				columns = line.rstrip().split("\t")

				if columns[0] != current_chr:
					print("     chr: %s"%(columns[0]))
					current_chr = columns[0]

				if columns[2] in set_esnp:
					dict_esnp[columns[2]] = [columns[0],int(columns[1])]

		stream_in.close()

		os.remove(vcf_in)

	# determine how many total eQTL there are
	total_eqtl = 0
	for esnp in df_qtl['variant_id'].tolist():
		if esnp in dict_esnp: total_eqtl += 1

	# 5 retrieve phenotype positions
	print("5. Retrieving ePhenotype positions...")
	set_epheno = set(df_qtl['phenotype_id'].tolist())

	stream_in = gzip.open(args.pheno, "r")

	dict_ephenotype = {}
	for line in stream_in:
		if isinstance(line, bytes) and not isinstance(line, str):
			line = line.decode()
		if line[0:1] != "#":
			columns = line.rstrip().split("\t")
			#Chr    start   end     ID
			if columns[3] in set_epheno:
				dict_ephenotype[columns[3]] = [columns[0],int(columns[1])+1]
	# print(dict_ephenotype)
	stream_in.close()

	# 6 calculate effect sizes
	print("6. Calculating eQTL effect sizes...")
	stream_vcf = open

	completed = 0
	t = time.time()

	stream_out = open(args.output, "w")
	headers = ['log2_aFC','log2_aFC_lower','log2_aFC_upper']
	if args.count_o == 1:
		headers += ['ref_allele_count','alt_allele_count']
	stream_out.write("\t".join(df_qtl.columns.tolist()+headers)+"\n")

	for index, row in df_qtl.iterrows():
		# now retrieve the genotypes for the snp
		# only for those individuals with phenotype data
		dict_geno = {}
		line_written = False
		allele_counts = None

		if row['variant_id'] in dict_esnp:
			if row['phenotype_id'] in dict_ephenotype:

				esnp = dict_esnp[row['variant_id']]
				records = tabix_vcf.fetch(esnp[0], esnp[1]-1, esnp[1])

				snp_found = 0
				for record in records:
					cols = record.rstrip().split("\t")
					if cols[2] == row['variant_id']:
						gt_index = cols[8].split(":").index(args.geno)
						snp_found = 1
						for sample in pheno_map.keys():
							sample_col = cols[vcf_map[sample]]
							dict_geno[sample] = sample_col.split(":")[gt_index]
				if snp_found == 0:
					print("	  WARNING: eSNP %s not found in VCF"%(row['variant_id']))
					continue

				# assume phenotype is within a megabase of SNP
				ephenotype = dict_ephenotype[row['phenotype_id']]
				# print(ephenotype)
				records = tabix_pheno.fetch(ephenotype[0], ephenotype[1]-1, ephenotype[1]+1)

				dict_pheno = {}

				for record in records:
					cols = record.rstrip().split("\t")
					if cols[3] == row['phenotype_id']:
						for sample in dict_geno.keys():
							if args.log_xform == 1 and args.log_base != 2:
								# if data has been log transformed but is not in base 2 convert it
								dict_pheno[sample] = float(cols[pheno_map[sample]]) * math.log(args.log_base,2)
							else:
								dict_pheno[sample] = float(cols[pheno_map[sample]])

				# make a dataframe with all covariates and genotype classes
				list_rows = []
				allele_counts = [0,0]

				for sample in dict_geno.keys():
					if args.geno == "GT":
						if "." not in dict_geno[sample]:	# only include samples w/ complete genotype data (no '.')
							list_rows.append([dict_geno[sample].count("1"),dict_pheno[sample]] + return_cov(sample))
							allele_counts[0] += dict_geno[sample].count("0")
							allele_counts[1] += dict_geno[sample].count("1")
					elif args.geno == "DS":
						list_rows.append([round(float(dict_geno[sample])),dict_pheno[sample]] + return_cov(sample))

				if len(list_rows) >= args.min_samps and min(allele_counts) >= args.min_alleles:		## Changed to only run effect size calc when more than minimum # samps w/ GT data and minimum number of observations for each allele
					df_test = pd.DataFrame(list_rows, columns=['geno','pheno']+["cov_"+x for x in df_cov[cov_id_col].tolist()])

					if args.matrix_o != None:
						df_test.to_csv(args.matrix_o+"/"+row['phenotype_id']+":"+row['variant_id']+".txt",sep="\t",index=False)

					# correct for covariates
					df_test = correct_covariates(df_test)

					esize = effect_size(df_test)
					line_out = row.tolist()+esize[0:3]
					if args.count_o == 1:
						line_out += allele_counts

					stream_out.write("\t".join(map(str,line_out))+"\n")
					line_written = True
				else:
					if len(list_rows) == 0:
						print("	  WARNING: no individual with genotype data for eQTL %s - %s"%(row['phenotype_id'],row['variant_id']))
					elif len(list_rows) < args.min_samps:
						print("	  WARNING: only %d individual(s) with genotype data for eQTL %s - %s"%(len(list_rows),row['phenotype_id'],row['variant_id']))
					elif min(allele_counts) < args.min_alleles:
						print("	  WARNING: only %d observations of minor allele for eQTL %s - %s"%(min(allele_counts),row['phenotype_id'],row['variant_id']))
			else:
				if row['variant_id'] != "nan" and args.chr == None:
					print("	  WARNING: positional information not found for ePhenotype %s"%(row['phenotype_id']))

			completed += 1

			if completed % 100 == 0:
				stream_out.flush()
				print("     COMPLETED %d of %d = %f in %d seconds"%(
				    completed, total_eqtl,
				    float(completed)/float(total_eqtl),
				    time.time()- t))
				t = time.time()
		else:
			if row['phenotype_id'] != "nan" and args.chr == None:
				print("	  WARNING: positional information not found for eSNP %s"%(row['variant_id']))
		
		# ensure that every eQTL has a line written so number of output lines = number of input eQTLs
		if ("sid_chr" in df_qtl.columns and (str(row['sid_chr']) == args.chr or args.chr == None)):
			if line_written == False:
				line_out = row.tolist()+[float('nan'),float('nan'),float('nan')]
				if args.count_o == 1:
					if allele_counts != None:
						line_out += allele_counts
					else:
						line_out += [float('nan'),float('nan')]
				stream_out.write("\t".join(map(str,line_out))+"\n")
	stream_out.close()

	duration = time.time() - start_time
	print("COMPLETED - total runtime was %d seconds"%(duration))


def return_cov(sample):
	global df_cov

	if sample in df_cov.columns:
		return(df_cov[sample].tolist())
	else:
		return([])


def sample_column_map(path, start_col=9, line_key="#CHR"):
	stream_in = gzip.open(path, "r")

	out_map = {}
	for line in stream_in:
		if isinstance(line, bytes) and not isinstance(line, str):
			line = line.decode()
		if line_key in line:
			line = line.rstrip().split("\t")
			for i in range(start_col,len(line)):
				out_map[line[i]] = i

			break

	stream_in.close()

	return(out_map)


def correct_covariates(df_test):
	global df_cov

	if len(df_cov.index) > 0:
		# correct for covariates
		# add genotype categorical covariates
		cov_homo_ref = [int(x == 0) for x in df_test['geno']]
		if sum(cov_homo_ref) > 0:
			df_test['cov_homo_ref'] = cov_homo_ref

		cov_homo_alt = [int(x == 2) for x in df_test['geno']]
		if sum(cov_homo_alt) > 0:
			df_test['cov_homo_alt'] = cov_homo_alt

		cov_ids = [x for x in df_test.columns if "cov_" in x]

		# convert categorical covariates to n-1 binary covariates
		new_cols = {}
		drop_cols = []

		for xcov in cov_ids:
			if df_test.dtypes[xcov] == object:
				values = list(set(df_test[xcov]))[1:]
				for xval in values:
					xname = xcov+"_"+xval
					new_cols[xname] = [int(x == xval) for x in df_test[xcov]]

				drop_cols.append(xcov)

		df_test.drop(drop_cols,axis=1,inplace=True)
		for xcov in new_cols.keys():
			df_test[xcov] = new_cols[xcov]
		cov_ids = [x for x in df_test.columns if "cov_" in x]

		# NOTE any variable that is a string will be treated as categorical - this is the same functionality as FASTQTL, so good
		# see: http://statsmodels.sourceforge.net/devel/example_formulas.html

		xformula = "pheno ~ "+"+".join(cov_ids)
		result = smf.ols(formula=xformula, data=df_test).fit()

		# use only significant (95% CI doesn't overlap 0) covariates to correct expression values
		# do not include intercept or genotypes in correction

		drop_covs = []
		for xcov in list(result.params.index):
			if xcov in df_test.columns:
				coefficient = result.params.loc[xcov]
				upper_ci = result.conf_int(0.05).loc[xcov][1]
				lower_ci = result.conf_int(0.05).loc[xcov][0]
				if (lower_ci <= 0 and upper_ci >= 0):
					drop_covs.append(xcov)

		# drop insignificant covariates
		df_test.drop(drop_covs, axis=1, inplace=True)
		cov_ids = [x for x in df_test.columns if "cov_" in x]

		# redo regression without insignificant covs
		if len(cov_ids) > 0:
			xformula = "pheno ~ "+"+".join(cov_ids)
			result = smf.ols(formula=xformula, data=df_test).fit()

			df_test_corrected = copy.deepcopy(df_test)
			for xcov in list(result.params.index):
				coefficient = result.params.loc[xcov]
				if xcov == "Intercept" or xcov == "cov_homo_ref" or xcov == "cov_homo_alt":
					df_test_corrected[xcov] = [0] * len(df_test_corrected.index)
				else:
					df_test_corrected[xcov] = [x * coefficient for x in df_test_corrected[xcov]]

			# add residual to dataframe
			df_test_corrected['pheno_cor'] = [row['pheno'] - sum(row[2:len(row)]) for index, row in df_test_corrected.iterrows()]

		else:
			# if none of the covariates are significant then just leave the values as is
			df_test_corrected = copy.deepcopy(df_test)
			df_test_corrected['pheno_cor'] = df_test_corrected['pheno']
	else:
		# covariates not provided
		df_test_corrected = copy.deepcopy(df_test)
		df_test_corrected['pheno_cor'] = df_test_corrected['pheno']


	return(df_test_corrected)


def effect_size(df_test):
	global args
	# calculate effect size
	esize = calculate_effect_size(df_test['geno'],df_test['pheno_cor'])

	# calculate 95% CI for effect size using BCa bootstrapping
	if args.boot > 0:
		try:
			ci = boot.ci((df_test['geno'].tolist(),df_test['pheno_cor'].tolist()), statfunction=calculate_effect_size, alpha=0.05, n_samples=args.boot, method="bca")
		except (IndexError,ValueError):		## ValueError added for calculating CI on one sample
			ci = [float('nan'),float('nan')]
	else:
		ci = [float('nan'),float('nan')]

	return([esize, ci[0],ci[1]])


def calculate_effect_size(genos,phenos):
	global args

	# in cases where there is only a single genotype in the data return nan
	if len(set(genos)) == 1:
		return(float('nan'))

	if args.log_xform == 1:
		# M5 - for log2 transformed data

		#1 need to prepare 4 estimates
		p_m = [
			np.mean(phenos[genos == 0]),
			np.mean(phenos[genos == 1]),
			np.mean(phenos[genos == 2])
		]

		log2ratio_M2M0 = bound_basic(p_m[2] - p_m[0], -args.ecap, args.ecap)
		log2ratio_M1M2 = bound_basic(p_m[1] - p_m[2], -0.9999999, args.ecap)
		log2ratio_M1M0 = bound_basic(p_m[1] - p_m[0], -1, args.ecap)

		p_delta = [
			float('nan'),
			math.pow(2,log2ratio_M2M0),
			float(1) / (math.pow(2,log2ratio_M1M2+1) - 1),
			math.pow(2,log2ratio_M1M0+1) - 1,
			None,
		]

		X = sm.add_constant(genos)
		result = sm.OLS(phenos,X).fit()
		result_coef = bound_basic(result.params[1]*2, -args.ecap, args.ecap)
		p_delta[4] = math.pow(2,result_coef)

		for x in range(1, 5):
			p_delta[x] = bound_basic(p_delta[x], math.pow(2,-args.ecap), math.pow(2,args.ecap))

		stdevs = {}

		# pick the estimate that minimizes residual variance
		for i in range(1,5):
			#stdevs[i] = numpy.std([yi - calculate_expected_expr(p_delta[i], xi) for xi, yi in zip(genos, phenos)])
			stdevs[i] = np.std(phenos - np.log2((2 - genos) + (p_delta[i] * genos)))

		min_delta = min([x for x in stdevs.values() if math.isnan(x) == False])
		use_delta = 0
		for delta in range(1,5):
			if stdevs[delta] == min_delta:
				use_delta = delta
				break
		p_delta[0] = float('nan')

		return(math.log(p_delta[use_delta],2))
	else:
		# linear regression on untransformed data
		X = sm.add_constant(genos)
		result = sm.OLS(phenos,X).fit()
		# ensure intercept is positive
		b0 = bound_basic(result.params[0], np.finfo(float).eps, float('inf'))
		# calculate the effect size
		use_delta = (float(2 * result.params[1]) / float(b0)) + 1

		# bound delta between caps
		if use_delta < 0:
			if result.params[1] > 0:
				use_delta = math.pow(2, args.ecap)
			else:
				use_delta = math.pow(2, -args.ecap)

		# bound effect size between -args.ecap and args.ecap in log space
		use_delta_log = math.log(use_delta, 2)
		use_delta_log_bounded = bound_basic(use_delta_log, -args.ecap, args.ecap)

		return(use_delta_log_bounded)


def bound_basic(x, l, h):
	y = min([x,h])
	y = max([y,l])
	if math.isnan(x) == True: 
		y = float('nan')
	return(y)


def calculate_expected_expr(delta, alt_alleles):
	if math.isnan(delta) == False:
		return(math.log((2 - alt_alleles) + (delta * alt_alleles), 2))
	else:
		return(float('nan'))


if __name__ == "__main__":
	main()
