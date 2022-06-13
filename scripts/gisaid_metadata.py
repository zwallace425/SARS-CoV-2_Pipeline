# Object used for initializing instance for GISAID metadata files. Takes the exact GISAID metadata file,
# in the format it's been since June 1st 2022, and transforms it into the dataframe format to analyze
# sequence prevalence or mutation prevalence dynamics with variant_scoring.py

# WARNING: GISAID metadata files are large.  We advise the user to store those files appropriately and to
# avoid frequent downloading of their data

import sys
import re
import pandas as pd
import pickle
import compress_json
from datetime import datetime
from variant_dynamics import VariantDynamics as vd

class GisaidMetadata(object):

	# Intialize GISAID metadata dataframe and perform QC on the data
	def __init__(self, gisaid_df, seq_length = 29400, n_content = 0.05, min_date = '2022-01-01'):
		columns = ['Virus name', 'Accession ID', 'Sequence length', 'AA Substitutions', 'Collection date', 'N-Content']
		gisaid_kept = gisaid_df[columns]
		gisaid_kept = gisaid_kept[(gisaid_kept['Sequence length'] > seq_length) & (gisaid_kept['N-Content'] < n_content)]
		gisaid_kept['Collection date'] = pd.to_datetime(gisaid_kept['Collection date'])
		first_current_month = datetime.today().replace(day=1).strftime('%Y-%m-%d')
		gisaid_kept = gisaid_kept[(gisaid_kept['Collection date'] < first_current_month) & (gisaid_kept['Collection date'] >= min_date)]
		gisaid_kept = gisaid_kept.drop_duplicates(subset = ['Accession ID'])
		gisaid_kept = gisaid_kept[(gisaid_kept['AA Substitutions'] != "()")]
		gisaid_kept = gisaid_kept.dropna()

		self.gisaid = gisaid_kept


	# The functions collects the covariate counts by region-date from the GISAID metadata file. If passed in
	# and exisiting list of accessions, it will only compute on the new accessions and update the exisiting covariate
	# prevalence JSON file as well as the region-date counts dictionary file
	"""Parameters:
		cov_accessions: A dictionary mapping GISAID acession to full covariate
		cov_prevalence: A dictionariy mapping covariate to region to date to counts
		cov_region_date_counts: A dictionary mapping region to date to total isolate counts
	"""
	def covariate_counts(self, cov_accessions = {}, cov_prevalence = {}, cov_region_date_counts = {}):

		data = self.gisaid

		if (cov_accessions != {} or cov_prevalence != {} or cov_region_date_counts != {}):
			if (cov_accessions == {} or cov_prevalence == {} or cov_region_date_counts == {}):
				sys.exit("Input Error: accessions, prevalence counts, and region-date counts must all be inputted or none inputed for GISAID analysis.")
		
		if cov_accessions != {}:
			cov_accessions = compress_json.load(cov_accessions)
			cov_prevalence = compress_json.load(cov_prevalence)
			cov_region_date_counts = compress_json.load(cov_region_date_counts)
			data = data[(not data['Acession ID'].isin(cov_accessions.keys()))]
	

		for ind in data.index:
			name = data['Virus name'][ind].split("/")
			region = name[1]
			acc = data['Accession ID'][ind]
			date = data['Collection date'][ind]
			date = date.strftime('%Y-%m-%d')

			if cov_region_date_counts.get(region):
				if cov_region_date_counts[region].get(date):
					cov_region_date_counts[region][date] += 1
				else:
					cov_region_date_counts[region][date] = 1
			else:
				cov_region_date_counts[region] = {}
				cov_region_date_counts[region][date] = 1

			cov = data['AA Substitutions'][ind]
			cov = cov[1:-1]
			cov = cov.split(",")
			cov.sort()
			cov = ",".join(cov)
			cov_accessions[acc] = cov

			if cov_prevalence.get(cov):
				if cov_prevalence[cov].get(region):
					if cov_prevalence[cov][region].get(date):
						cov_prevalence[cov][region][date] += 1
					else:
						cov_prevalence[cov][region][date] = 1
				else:
					cov_prevalence[cov][region] = {}
					cov_prevalence[cov][region][date] = 1
			else:
				cov_prevalence[cov] = {}
				cov_prevalence[cov][region] = {}
				cov_prevalence[cov][region][date] = 1

		compress_json.dump(cov_accessions, "data/gisaid_cov_accessions.json.gz")
		compress_json.dump(cov_prevalence, "data/gisaid_cov_prevalence.json.gz")
		compress_json.dump(cov_region_date_counts, "data/gisaid_cov_region_date_counts.json.gz")

		return (vd.variant_dict_to_df(cov_prevalence, 'gisaid_cov'), 
        	vd.region_date_dict_to_df(cov_region_date_counts, 'gisaid_cov'))
	

    # The functions collects the AA mutation counts by region-date from the GISAID metadata file. If passed in
	# and exisiting list of accessions, it will only compute on the new accessions and update the exisiting AA mutations
	# prevalence JSON file as well as the region-date counts dictionary file
	"""Parameters:
		aa_accessions: A dictionary mapping GISAID acession to full covariate (like in the function above but this is only supposed
			to be for accessions that we've computed AA mutation counts)
		aa_prevalence: A dictionariy mapping AA mutation to region to date to counts
		aa_region_date_counts: A dictionary mapping region to date to total isolate counts
	"""
	def mutation_counts(self, aa_accessions = {}, aa_prevalence = {}, aa_region_date_counts = {}):

		data = self.gisaid

		if (aa_accessions != {} or aa_prevalence != {} or aa_region_date_counts != {}):
			if (aa_accessions == {} or aa_prevalence == {} or aa_region_date_counts == {}):
				sys.exit("Input Error: accessions, prevalence counts, and region-date counts must all be inputted or none inputed for GISAID analysis.")
		
		if aa_accessions != {}:
			aa_accessions = compress_json.load(aa_accessions)
			aa_prevalence = compress_json.load(aa_prevalence)
			aa_region_date_counts = compress_json.load(aa_region_date_counts)
			data = data[(not data['Acession ID'].isin(aa_accessions.keys()))]


		for ind in data.index:
			name = data['Virus name'][ind].split("/")
			region = name[1]
			acc = data['Accession ID'][ind]
			date = data['Collection date'][ind]
			date = date.strftime('%Y-%m-%d')

			if aa_region_date_counts.get(region):
				if aa_region_date_counts[region].get(date):
					aa_region_date_counts[region][date] += 1
				else:
					aa_region_date_counts[region][date] = 1
			else:
				aa_region_date_counts[region] = {}
				aa_region_date_counts[region][date] = 1

			cov = data['AA Substitutions'][ind]
			cov = cov[1:-1]
			aa_accessions[acc] = cov
			muts = cov.split(",")
			aa_prevalence = vd.aa_prevalence_main(aa_prevalence, muts, region, date)

		compress_json.dump(aa_accessions, "data/gisaid_aa_accessions.json.gz")
		compress_json.dump(aa_prevalence, "data/gisaid_aa_prevalence.json.gz")
		compress_json.dump(aa_region_date_counts, "data/gisaid_aa_region_date_counts.json.gz")

		return (vd.variant_dict_to_df(aa_prevalence, 'gisaid_aa'), 
			vd.region_date_dict_to_df(aa_region_date_counts, 'gisaid_aa'))

	

	# The functions collects the Spike covariate counts by region-date from the GISAID metadata file. If passed in
	# and exisiting list of accessions, it will only compute on the new accessions and update the exisiting Spike covariate
	# prevalence JSON file as well as the region-date counts dictionary file
	"""Parameters:
		spike_accessions: A dictionary mapping GISAID acession to full covariate (like in the function above but this is only supposed
			to be for accessions that we've computed Spike covariate counts)
		spike_prevalence: A dictionariy mapping Spike covariate to region to date to counts
		spike_region_date_counts: A dictionary mapping region to date to total isolate counts
	"""
	def spike_counts(self, spike_accessions = {}, spike_prevalence = {}, spike_region_date_counts = {}):

		def num_sort(text):
			return list(map(int, re.findall(r'\d+', text)))[0]

		data = self.gisaid

		if (spike_accessions != {} or spike_prevalence != {} or spike_region_date_counts != {}):
			if (spike_accessions == {} or spike_prevalence == {} or spike_region_date_counts == {}):
				sys.exit("Input Error: accessions, prevalence counts, and region-date counts must all be inputted or none inputed for GISAID analysis.")
		
		if spike_accessions != {}:
			spike_accessions = compress_json.load(spike_accessions)
			spike_prevalence = compress_json.load(spike_prevalence)
			spike_region_date_counts = compress_json.load(spike_region_date_counts)
			data = data[(not data['Acession ID'].isin(spike_accessions.keys()))]

		for ind in data.index:
			name = data['Virus name'][ind].split("/")
			region = name[1]
			acc = data['Accession ID'][ind]
			date = data['Collection date'][ind]
			date = date.strftime('%Y-%m-%d')

			if spike_region_date_counts.get(region):
				if spike_region_date_counts[region].get(date):
					spike_region_date_counts[region][date] += 1
				else:
					spike_region_date_counts[region][date] = 1
			else:
				spike_region_date_counts[region] = {}
				spike_region_date_counts[region][date] = 1

			cov = data['AA Substitutions'][ind]
			cov = cov[1:-1]
			spike_accessions[acc] = cov
			muts = cov.split(",")
			
			spike_cov = []
			for mut in muts:
				if "Spike" in mut:
					spike_cov.append(mut)
			spike_cov.sort(key=num_sort)
			spike_cov = ",".join(spike_cov)

			if spike_prevalence.get(spike_cov):
				if spike_prevalence[spike_cov].get(region):
					if spike_prevalence[spike_cov][region].get(date):
						spike_prevalence[spike_cov][region][date] += 1
					else:
						spike_prevalence[spike_cov][region][date] = 1
				else:
					spike_prevalence[spike_cov][region] = {}
					spike_prevalence[spike_cov][region][date] = 1
			else:
				spike_prevalence[spike_cov] = {}
				spike_prevalence[spike_cov][region] = {}
				spike_prevalence[spike_cov][region][date] = 1

		compress_json.dump(spike_accessions, "data/gisaid_spike_accessions.json.gz")
		compress_json.dump(spike_prevalence, "data/gisaid_spike_prevalence.json.gz")
		compress_json.dump(spike_region_date_counts, "data/gisaid_spike_region_date_counts.json.gz")

		return (vd.variant_dict_to_df(spike_prevalence, 'gisaid_spike'), 
			vd.region_date_dict_to_df(spike_region_date_counts, 'gisaid_spike'))














