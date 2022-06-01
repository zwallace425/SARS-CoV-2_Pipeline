# Object for analyzing SARS-CoV-2 variants based on epedimiological data and prior experimental knowledge

import sys
import pickle
import numpy as np
import pandas as pd
import compress_json
from datetime import datetime
import variant_dynamics as vd

class VariantAnalysis(object):

	# This function will analyze the saved compressed JSON that was created in spatial_temporal_prevalence in
	# variant_dynamics.py as well as uncompress the region_date json file
	@staticmethod
	def analyze_from_saved_dict(prevalence_dict_json, region_date_dict_json, interval = 'W', world = False):
			
		prevalence_dict = compress_json.load(prevalence_dict_json)
		region_date_dict = compress_json.load(region_date_dict_json)

		prevalence_df = []
		for variant in prevalence_dict.keys():
			for region in prevalence_dict[variant].keys():
				for date in prevalence_dict[variant][region].keys():
					prev = prevalence_dict[variant][region][date]
					dto = datetime.strptime(date, '%Y_%m_%d').date()
					data = pd.DataFrame({'Variant': [variant], 'Region': [region], 'Date': [dto], 'Count': [prev]})
					prevalence_df.append(data)

		region_date_df = []
		for region in region_date_dict.keys():
			for date in region_date_dict[region].keys():
				count = region_date_dict[region][date]
				dto = datetime.strptime(date, '%Y_%m_%d').date()
				data = pd.DataFrame({'Region': [region], 'Date': [dto], 'Count': [count]})
				region_date_df.append(data)

		prevalence_df = pd.concat(prevalence_df, ignore_index = True)
		region_date_df = pd.concat(region_date_df, ignore_index = True)

		return(VariantAnalysis.analyze_dynamics(prevalence_df, region_date_df, interval, world))


	# This function will analyze the saved pickled dataframe that was created in dictionary_to_dataframe
	# in variant_dynamics.py and also open the region_date pickeled dataframe
	@staticmethod
	def analyze_from_saved_df(prevalence_df_pkl, region_date_df_pkl, interval = 'W', world = False):
		with open(prevalence_df_pkl, 'rb') as handle:
			prevalence_df = pickle.load(handle)

		with open(region_date_df_pkl, 'rb') as handle:
			region_date_df = pickle.load(handle)

		return(VariantAnalysis.analyze_dynamics(prevalence_df, region_date_df, interval, world))


	# This will manipulate the dataframe generated from analyze_covariates_from_output (or likewise for the AA function)
	# that has the covariate/aa mutation counts by date and region into a dataframe grouped by a specified time interval
	# by region, or globaly if "world" is True.  So if a user wants to analyze variant dynamics by week, this will group
	# all variant counts into its appropriate week bin by region and then use those aggregations to compute the prevalence
	# ratio, growth, and jerk by regions by date.  If "world" is True that the dataframe will represent dynamics of each
	# variant globally by date. 
	@staticmethod
	def analyze_dynamics(prevalence_df, region_date_df, interval = 'W', world = False):

		prevalence_df['Date'] = pd.to_datetime(prevalence_df['Date'])
		region_date_df['Date'] = pd.to_datetime(region_date_df['Date'])
		date_group = pd.Grouper(key = 'Date', freq = interval)

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


		return(prevalence_df)
			


if __name__ == "__main__":

	reference = sys.argv[1]
	sequences = sys.argv[2]
	interval = sys.argv[3]
	world = sys.argv[4]


	if (sys.argv[4] == "T"):
		world = True
	elif (sys.argv[4] == "F"):
		world = False

	cov_prev, aa_prev, region_date_counts = vd.VariantDynamics.spatial_temporal_prevalence(reference, sequences)

	cov_prev_df = VariantAnalysis.analyze_dynamics(cov_prev, region_date_counts, interval, world)
	aa_prev_df = VariantAnalysis.analyze_dynamics(aa_prev, region_date_counts, interval, world)

	print(cov_prev_df)
	print('\n')
	print(aa_prev_df)

