# Object for provding an analyzable dataframe of epedimologcial dynamics for SARS-Cov-2 variants.
# The general return from this object is a dataframe of variant prevalence, growth, and jerk by 
# date bin (either a single day, week, 2-weeks, or month) per country.  It requires the dataframe
# from VariantDynamics

import sys
import pickle
import numpy as np
import pandas as pd
import compress_json
from datetime import datetime
import variant_dynamics as vd

class VariantAnalysis(object):

	# This function will analyze the saved compressed JSON that was created in spatial_temporal_prevalence in
	# variant_dynamics.py as well as uncompress the region_date json .
	"""Parameters:
		prevalence_df: Dictionary of aa mutation or covariate counts by region and date stored as compressed JSON
		region_date_df: Dictionary of region-date isolate counts stored as compressed JSON
		period: Period of time by which to aggregate counts by date and compute variant prevalence from period to period.
			This paramter must either be 'D' (day), 'W' (week), '2W' (bi-weekly), or 'M' (month)
		world: Boolean to determine whether to aggregate prevalence globally or leave by region 
	"""
	@staticmethod
	def analyze_from_saved_dict(prevalence_dict_json, region_date_dict_json, period = '2W', world = False):
			
		prevalence_dict = compress_json.load(prevalence_dict_json)
		region_date_dict = compress_json.load(region_date_dict_json)

		prevalence_df = []
		for variant in prevalence_dict.keys():
			for region in prevalence_dict[variant].keys():
				for date in prevalence_dict[variant][region].keys():
					prev = prevalence_dict[variant][region][date]
					dto = datetime.strptime(date, '%Y-%m-%d').date()
					data = pd.DataFrame({'Variant': [variant], 'Region': [region], 'Date': [dto], 'Count': [prev]})
					prevalence_df.append(data)

		region_date_df = []
		for region in region_date_dict.keys():
			for date in region_date_dict[region].keys():
				count = region_date_dict[region][date]
				dto = datetime.strptime(date, '%Y-%m-%d').date()
				data = pd.DataFrame({'Region': [region], 'Date': [dto], 'Count': [count]})
				region_date_df.append(data)

		prevalence_df = pd.concat(prevalence_df, ignore_index = True)
		region_date_df = pd.concat(region_date_df, ignore_index = True)

		return(VariantAnalysis.analyze_dynamics(prevalence_df, region_date_df, period, world))


	# This function will analyze the saved pickled dataframe that was created in dictionary_to_dataframe
	# in variant_dynamics.py and also open the region_date pickeled dataframe
	"""Parameters:
		prevalence_df: Pickeled dataframe of aa mutation or covariate counts by region and date
		region_date_df: Pickeled dataframe of region-date isolate counts
		period: Period of time by which to aggregate counts by date and compute variant prevalence from period to period.
			This paramter must either be 'D' (day), 'W' (week), '2W' (bi-weekly), or 'M' (month)
		world: Boolean to determine whether to aggregate prevalence globally or leave by region 
	"""
	@staticmethod
	def analyze_from_saved_df(prevalence_df_pkl, region_date_df_pkl, period = '2W', world = False):
		
		with open(prevalence_df_pkl, 'rb') as handle:
			prevalence_df = pickle.load(handle)

		with open(region_date_df_pkl, 'rb') as handle:
			region_date_df = pickle.load(handle)

		return(VariantAnalysis.analyze_dynamics(prevalence_df, region_date_df, period, world))


	# This function will manipulate the protein covariates prevalence dictionary, which has the date-region prevalence
	# counts of each protein specific covariates generated from spatial_temporal_prevalence in VariantDynamics and return
	# this as a protein covariate specific dataframe to be analyzeed in analyze dynamics.  This function is called when 
	# a user wants to analyze dynamics of just spike covariates, or any othe protein, instead of the full variant constellation.
	"""Parameters:
		prot_cov_prev_dict: Dictionary of protein-specific covariate region-date counts
		region_date_df: Dataframe of region-date isolate counts
		period: Period of time by which to aggregate counts by date and compute variant prevalence from period to period.
			This paramter must either be 'D' (day), 'W' (week), '2W' (bi-weekly), or 'M' (month)
		world: Boolean to determine whether to aggregate prevalence globally or leave by region 
	"""
	@staticmethod
	def analyze_protein_covariates(prot_cov_prev_dict, region_date_df, protein, period = '2W', world = False):
		
		if (prot_cov_prev_dict.get(protein) == None):
			sys.exit("Invalid SARS-CoV-2 protein.  Please enter a valid non-structural protein name of accessory protein name.")

		prevalence_dict = prot_cov_prev_dict[protein]

		prevalence_df = []
		for variant in prevalence_dict.keys():
			for region in prevalence_dict[variant].keys():
				for date in prevalence_dict[variant][region].keys():
					prev = prevalence_dict[variant][region][date]
					dto = datetime.strptime(date, '%Y-%m-%d').date()
					data = pd.DataFrame({'Variant': [variant], 'Region': [region], 'Date': [dto], 'Count': [prev]})
					prevalence_df.append(data)

		prevalence_df = pd.concat(prevalence_df, ignore_index = True)

		return(VariantAnalysis.analyze_dynamics(prevalence_df, region_date_df, period, world))


	# This will manipulate the dataframe generated from spatial_temporal_prevalence in VariantDynamics, cov_prevalence or aa_prevalence,
	# that has the covariate/aa mutation counts by date and region into a dataframe grouped by a specified time period
	# by region, or globaly if "world" is True.  So if a user wants to analyze variant dynamics by week, this will group
	# all variant counts into its appropriate week bin by region and then use those aggregations to compute the prevalence
	# ratio, growth, and jerk by regions by date.  If "world" is True that the dataframe will represent dynamics of each
	# variant globally by date.
	"""Parameters:
		prevalence_df: Dataframe of aa mutation or covariate counts by region and date
		region_date_df: Dataframe of region-date isolate counts
		period: Period of tme by which to aggregate counts by date and compute variant prevalence from period to period.
			This paramter must either be 'D' (day), 'W' (week), '2W' (bi-weekly), or 'M' (month)
		world: Boolean to determine whether to aggregate prevalence globally or leave by region 
	"""
	@staticmethod
	def analyze_dynamics(prevalence_df, region_date_df, period = '2W', world = False):

		print("Aggregating counts and computing growth dynamics over time per region ...")

		prevalence_df['Date'] = pd.to_datetime(prevalence_df['Date'])
		region_date_df['Date'] = pd.to_datetime(region_date_df['Date'])
		date_group = pd.Grouper(key = 'Date', freq = period)

		if world:
			time_counts = region_date_df.groupby([date_group])['Count'].sum().reset_index().sort_values('Date')
			region_counts = prevalence_df.groupby(['Variant', 'Region', date_group]).size().groupby(['Variant', 'Date']).size().to_frame(
				name = 'Total Regions').reset_index().sort_values('Date')
			prevalence_df = prevalence_df.groupby(['Variant', date_group])['Count'].sum().reset_index().sort_values('Date').merge(
				time_counts, on = ['Date']).rename(columns = {'Count_x': 'Variant Count','Count_y': 'Isolates Count'})
			prevalence_df['Prevalence'] = (prevalence_df['Variant Count']/prevalence_df['Isolates Count'])
			prevalence_df['Growth'] = prevalence_df.groupby(['Variant'])['Prevalence'].pct_change().add(1)
			prevalence_df['Jerk'] = prevalence_df.groupby(['Variant'])['Growth'].diff()
			prevalence_df = prevalence_df.fillna(0).merge(region_counts, on = ['Variant', 'Date'])
		else:
			region_time_count = region_date_df.groupby(['Region', date_group])['Count'].sum().reset_index().sort_values('Date')
			prevalence_df = prevalence_df.groupby(['Variant', 'Region', date_group])['Count'].sum().reset_index().sort_values(['Date']).merge(
				region_time_count, on = ['Region', 'Date']).rename(columns = {'Count_x': 'Variant Count', 'Count_y': 'Isolates Count'})
			prevalence_df['Prevalence'] = (prevalence_df['Variant Count']/prevalence_df['Isolates Count'])
			prevalence_df['Growth'] = prevalence_df.groupby(['Variant', 'Region'])['Prevalence'].pct_change().add(1)
			prevalence_df['Jerk'] = prevalence_df.groupby(['Variant', 'Region'])['Growth'].diff()
			prevalence_df = prevalence_df.fillna(0)

		prevalence_df['Variant'].replace('', np.nan, inplace=True)
		prevalence_df.dropna(subset = ['Variant'], inplace = True)

		print("Done computing growth dynamics")
		print('\n')

		return(prevalence_df)
			


if __name__ == "__main__":

	reference = sys.argv[1]
	sequences = sys.argv[2]
	period = sys.argv[3]
	world = sys.argv[4]
	protein = sys.argv[5]


	if (sys.argv[4] == "T"):
		world = True
	elif (sys.argv[4] == "F"):
		world = False

	cov_prev, aa_prev, prot_cov_prev, region_date_counts = vd.VariantDynamics.spatial_temporal_prevalence(reference, sequences)

	cov_prev_df = VariantAnalysis.analyze_dynamics(cov_prev, region_date_counts, period, world)
	aa_prev_df = VariantAnalysis.analyze_dynamics(aa_prev, region_date_counts, period, world)
	prot_cov_pref_df = VariantAnalysis.analyze_protein_covariates(prot_cov_prev, region_date_counts, protein, period, world)

	print(cov_prev_df)
	print('\n')
	print(aa_prev_df)
	print('\n')
	print(prot_cov_pref_df)

