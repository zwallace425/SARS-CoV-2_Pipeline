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
	def __init__(self, gisaid_df, seq_length = 29400, n_content = 0.01, min_date = '2022-01-01', max_date = None, WHO = None, PANGO = None):
		
		print("Preparing GISAID Metadata ...")
		
		if WHO and PANGO:
			sys.exit("Input Error: WHO name and PANGO Lineage name not allowed to be iputted together.")
		
		columns = ['Virus name', 'Accession ID', 'Sequence length', 'AA Substitutions', 'Variant', 'Pango lineage','Collection date', 'N-Content']
		gisaid_kept = gisaid_df[columns]
		gisaid_kept = gisaid_kept[(gisaid_kept['Sequence length'] > seq_length) & (gisaid_kept['N-Content'] < n_content)]
		gisaid_kept['Collection date'] = pd.to_datetime(gisaid_kept['Collection date'])
		if max_date:
			max_month = datetime.strptime(max_date, '%Y-%m-%d')
			max_month = max_month.replace(day=1).strftime('%Y-%m-%d')
		else:
			max_month = datetime.today().replace(day=1).strftime('%Y-%m-%d')
		gisaid_kept = gisaid_kept[(gisaid_kept['Collection date'] < max_month) & (gisaid_kept['Collection date'] >= min_date)]
		gisaid_kept = gisaid_kept.drop_duplicates(subset = ['Accession ID'])
		gisaid_kept = gisaid_kept[(gisaid_kept['AA Substitutions'] != "()")]
		gisaid_kept = gisaid_kept.dropna()

		print("Done preparing GISAID data")
		print('\n')

		self.gisaid = gisaid_kept
		self.who = WHO
		self.pango = PANGO


	# The function collects the covariate counts by region-date from the GISAID metadata file.
	def covariate_counts(self):

		print("Aggregating covariate region-date counts ...")

		data = self.gisaid
		cov_prevalence = {}
		cov_region_date_counts = {}
		cov_pango = {}

		for ind in data.index:
			name = data['Virus name'][ind].split("/")
			region = name[1]
			date = data['Collection date'][ind]
			date = date.strftime('%Y-%m-%d')
			pango_lineage = data['Pango lineage'][ind]

			if cov_region_date_counts.get(region):
				if cov_region_date_counts[region].get(date):
					cov_region_date_counts[region][date] += 1
				else:
					cov_region_date_counts[region][date] = 1
			else:
				cov_region_date_counts[region] = {}
				cov_region_date_counts[region][date] = 1

			if self.who:
				data_type = 'gisaid_'+self.who+'_cov'
				clade = data['Variant'][ind]
				if self.who in clade:
					cov = data['AA Substitutions'][ind]
				else:
					continue
			elif self.pango:
				data_type = 'gisaid_'+self.pango+'_cov'
				if self.pango in pango_lineage:
					cov = data['AA Substitutions'][ind]
				else:
					continue
			else:
				data_type = 'gisaid_cov'
				cov = data['AA Substitutions'][ind]

			cov = cov[1:-1]
			cov = cov.split(",")
			cov.sort()
			cov = ",".join(cov)
			cov_pango[cov] = pango_lineage

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

		print("Done aquiring region-date counts")
		print('\n')

		return (vd.variant_dict_to_df(cov_prevalence, data_type), 
        	vd.region_date_dict_to_df(cov_region_date_counts, 'gisaid'),
        	GisaidMetadata.cov_pango_to_df(cov_pango))
	

	# The functions collects the AA mutation counts by region-date from the GISAID metadata file.
	"""Parameters:
		protein: Optional protein input to only compute mutations for
	"""
	def mutation_counts(self, protein = ""):

		print("Aggregating mutation region-date counts ...")

		data = self.gisaid
		aa_prevalence = {}
		aa_region_date_counts = {}

		for ind in data.index:
			name = data['Virus name'][ind].split("/")
			region = name[1]
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

			if self.who:
				data_type = 'gisaid_'+self.who+'_aa'
				clade = data['Variant'][ind]
				if self.who in clade:
					cov = data['AA Substitutions'][ind]
				else:
					continue
			elif self.pango:
				data_type = 'gisaid_'+self.pango+'_aa'
				pango_lineage = data['Pango lineage'][ind]
				if self.pango in pango_lineage:
					cov = data['AA Substitutions'][ind]
				else:
					continue
			else:
				data_type = 'gisaid_aa'
				cov = data['AA Substitutions'][ind]

			cov = cov[1:-1]
			muts = cov.split(",")
			muts = [ m for m in muts if protein+'_' in m ]
			aa_prevalence = vd.aa_prevalence_main(aa_prevalence, muts, region, date)

		print("Done aquiring region-date counts")
		print('\n')

		return (vd.variant_dict_to_df(aa_prevalence, data_type), 
			vd.region_date_dict_to_df(aa_region_date_counts, 'gisaid'))



	# The functions collects the covariate counts by region-date from the GISAID metadata file for a specific protein.
	"""Parameters:
		protein: Inputted protein to compute the covariates for
	"""
	def protein_counts(self, protein):

		print("Aggregating protein-specific covariate region-date counts ...")

		def num_sort(text):
			return list(map(int, re.findall(r'\d+', text)))[0]

		data = self.gisaid
		prot_prevalence = {}
		prot_region_date_counts = {}

		for ind in data.index:
			name = data['Virus name'][ind].split("/")
			region = name[1]
			date = data['Collection date'][ind]
			date = date.strftime('%Y-%m-%d')

			if prot_region_date_counts.get(region):
				if prot_region_date_counts[region].get(date):
					prot_region_date_counts[region][date] += 1
				else:
					prot_region_date_counts[region][date] = 1
			else:
				prot_region_date_counts[region] = {}
				prot_region_date_counts[region][date] = 1

			if self.who:
				data_type = 'gisaid_'+self.who+'_'+protein+'_cov'
				clade = data['Variant'][ind]
				if self.who in clade:
					cov = data['AA Substitutions'][ind]
				else:
					continue
			elif self.pango:
				data_type = 'gisaid_'+self.pango+'_'+protein+'_cov'
				pango_lineage = data['Pango lineage'][ind]
				if self.pango in pango_lineage:
					cov = data['AA Substitutions'][ind]
				else:
					continue
			else:
				data_type = 'gisaid_'+protein+'_cov'
				cov = data['AA Substitutions'][ind]

			cov = cov[1:-1]
			muts = cov.split(",")
			
			prot_cov = []
			for mut in muts:
				m = mut.split("_")[0]
				if protein == m:
					prot_cov.append(mut)
			prot_cov.sort(key=num_sort)
			prot_cov = ",".join(prot_cov)

			if prot_prevalence.get(prot_cov):
				if prot_prevalence[prot_cov].get(region):
					if prot_prevalence[prot_cov][region].get(date):
						prot_prevalence[prot_cov][region][date] += 1
					else:
						prot_prevalence[prot_cov][region][date] = 1
				else:
					prot_prevalence[prot_cov][region] = {}
					prot_prevalence[prot_cov][region][date] = 1
			else:
				prot_prevalence[prot_cov] = {}
				prot_prevalence[prot_cov][region] = {}
				prot_prevalence[prot_cov][region][date] = 1

		print("Done acquiring protein-specific region-date counts")
		print('\n')

		return (vd.variant_dict_to_df(prot_prevalence, data_type),
			vd.region_date_dict_to_df(prot_region_date_counts, 'gisaid'))


	# This function collects PANGO Lineage counts by region and date from the GISAID Metadata file.
	# NOTE: The data returned from this function is not used in any of our scoring heuristics, as the
	# Emerging Lineage Score actually requires covariate counts while mapping covariates to the lineage.
	# However, this function is needed to collect the data for plotting lineage growth over time, as
	# computed in variant_plots.py
	def lineage_counts(self):

		print("Aggregating PANGO Lineage region-date counts ...")
		
		data = self.gisaid
		lineage_prevalence = {}
		lineage_region_date_counts = {}

		for ind in data.index:
			lineage = data['Pango lineage'][ind]
			name = data['Virus name'][ind].split("/")
			region = name[1]
			date = data['Collection date'][ind]
			date = date.strftime('%Y-%m-%d')

			if lineage_region_date_counts.get(region):
				if lineage_region_date_counts[region].get(date):
					lineage_region_date_counts[region][date] += 1
				else:
					lineage_region_date_counts[region][date] = 1
			else:
				lineage_region_date_counts[region] = {}
				lineage_region_date_counts[region][date] = 1

			if lineage_prevalence.get(lineage):
				if lineage_prevalence[lineage].get(region):
					if lineage_prevalence[lineage][region].get(date):
						lineage_prevalence[lineage][region][date] += 1
					else:
						lineage_prevalence[lineage][region][date] = 1
				else:
					lineage_prevalence[lineage][region] = {}
					lineage_prevalence[lineage][region][date] = 1
			else:
				lineage_prevalence[lineage] = {}
				lineage_prevalence[lineage][region] = {}
				lineage_prevalence[lineage][region][date] = 1

		print("Done acquiring PANGO Lineage region-date counts")
		print('\n')

		return (vd.variant_dict_to_df(lineage_prevalence, 'lineage'), 
			vd.region_date_dict_to_df(lineage_region_date_counts, 'gisaid'))


	# This function is just a helper for coverting the mapping of full variant constellations to pango lineage into a dataframe
	@staticmethod
	def cov_pango_to_df(cov_pango_dict):

		cov_pango_df = []
		for cov in cov_pango_dict.keys():
			pango = cov_pango_dict[cov]
			data = pd.DataFrame({'PANGO Lineage': [pango], 'Variant': [cov]})
			cov_pango_df.append(data)

		cov_pango_df = pd.concat(cov_pango_df, ignore_index = True)

		return cov_pango_df












