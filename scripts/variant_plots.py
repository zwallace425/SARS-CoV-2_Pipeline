# Object for computing growth dynamics plot of SARS-CoV-2 variants. This only accepts a dataframe in the format returned
# from analyze_dynamics in the VariantAnalysis object (same format required for VariantScoring), which is ...

"""	Variant 		Region		Date 		Variant Count 		Isolates Count 		Prevalence 		Growth 		Jerk
	X				X			YYYY-MM-DD	X 					X 					X 				X 			X
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class VariantPlots(object):

	# Intialize the variant dynamics dataframe to only contain the last six timestamps of data
	def __init__(self, variants_df, interval = 6, region = 'World', period = 'M'):

		dates = variants_df.sort_values(by = 'Date', ascending = False).drop_duplicates(subset = ['Date']).head(interval)
		dates = list(dates['Date'])
		dynamics_df = variants_df[(variants_df['Date'].isin(dates))]

		if region != 'World':
			try:
				dynamics_df = dynamics_df[(dynamics_df['Region'] == region)]
			except:
				raise Exception("Invalid dynamics dataframe. No column 'Region'.")
		
		self.dynamics_df = dynamics_df
		self.dates = dates
		self.region = region
		self.period = period


	# This function plots the growth of PANGO Lineages over the course of six months in a proportional stack plot. If no list 
	# of lineages is supplied, then it plots the 9 lineages with the highest prevalence ratios over the last six months, with a tenth
	# being all others. This function assumes that the initialized instance was a dataframe of lineage dynamic either globally or 
	# by region.
	"""Parameters:
		lineages: Optional list of lineages to specifically plot
	"""
	def plot_lineages(self, lineages = []):

		data = self.dynamics_df

		if lineages:
			try:
				data = data[(data['Variant'].isin(lineages))]
			except:
				raise Exception("Invalid PANGO Lineages input. Could not find all user supplied lineages.")

		top_lineages = data.sort_values(by = "Prevalence", ascending = False).drop_duplicates(subset = "Variant").head(14)
		top_lineages = list(top_lineages['Variant'])
		data = data[data['Variant'].isin(top_lineages)]
		plot_data = pd.DataFrame({'Lineage': top_lineages})
		for date in self.dates:
			date_data = data[data['Date'] == date]
			year_month = date.to_period('M').strftime("%Y-%m")
			date_growth_data = pd.DataFrame({'Lineage': list(date_data['Variant']), year_month: list(date_data['Prevalence'])})
			plot_data = plot_data.join(date_growth_data.set_index('Lineage'), on='Lineage')

		plot_data = plot_data.fillna(0)
		plot_data = plot_data.set_index('Lineage')
		plot_data.loc['Other'] = 1 - plot_data.sum(axis = 0)
		plot_data = plot_data[plot_data.columns[::-1]]
		x_val = np.asarray(plot_data.columns.tolist())
		indices = plot_data.index.tolist()
		plt.stackplot(x_val, plot_data.loc[indices].values, labels = indices)
		plt.suptitle('A Time Course of Lineage Proportions from Sequences Isolated in '+self.region, fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Sequence Prevalence Ratio (Variant Count/Isolates Count)')
		plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.10), ncol = 5, fancybox = True, shadow = True)
		plt.show()


	# This function plots the growth of covariates over the course of six months in a proportional stack plots. If not list of covariates
	# is supplied, then it plots the 10 covariates with the highest prevalence ratio.  The __init__ function took in a pre-established 
	# variant dynamics dataframe that could be protein specific covariates, PANGO/WHO specific covariates, or a combination of
	# protein-pango, protein-who covariates; hence, to accurately label the visualization, the protein name or PANGO/WHO clade will need
	# to be passed in to the function.
	"""Parameters:
		covariates: Optional list mutation to specificially plot
		protein: Optional protein name, but recomended if initialized dynamics dataframe is protein specific to  label plot
		PANGO: Optional Pango lineage name, but recomended if initialized dynamics dataframe is lineage specific to label plot
		WHO: Optional WHO name, but recommended if initialized dynamics dataframe is WHO name specific to label plot
	"""
	def plot_covariates(self, protein = "", PANGO = "", WHO = "", covariates = []):

		data = self.dynamics_df

		if WHO and PANGO:
			sys.exit("Input Error: WHO name and PANGO Lineage name not allowed to be iputted together.")

		if covariates:
			try:
				data = data[(data['Variant'].isin(covariates))]
			except:
				raise Exception("Invalid covariates input. Could not find all user supplied covariates.")

		if protein:
			data['Variant'] = data['Variant'].str.replace(protein+"_", "")

		top_covariates = data.sort_values(by = "Prevalence", ascending = False).drop_duplicates(subset = "Variant").head(15)
		top_covariates = list(top_covariates['Variant'])
		data = data[data['Variant'].isin(top_covariates)]
		plot_data = pd.DataFrame({'Variant': top_covariates})
		for date in self.dates:
			date_data = data[data['Date'] == date]
			year_month = date.to_period('M').strftime("%Y-%m")
			date_growth_data = pd.DataFrame({'Variant': list(date_data['Variant']), year_month: list(date_data['Prevalence'])})
			plot_data = plot_data.join(date_growth_data.set_index('Variant'), on='Variant')

		plot_data = plot_data.fillna(0)
		plot_data = plot_data.set_index('Variant')
		plot_data = plot_data[plot_data.columns[::-1]]
		x_val = np.asarray(plot_data.columns.tolist())
		indices = plot_data.index.tolist()
		plt.stackplot(x_val, plot_data.loc[indices].values, labels = indices)
		plt.suptitle('A Time Course of '+PANGO+WHO+' '+protein+' Covariate Prevalence Ratios from Sequences Isolated in '+self.region, fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Sequence Prevalence Ratio (Variant Count/Isolates Count)')
		plt.legend(loc = 'upper center', prop = {'size': 7}, bbox_to_anchor = (0.5, 1.10))
		plt.show()


	# This function plots the growth of amino acid mutations over the course of six months in a single plot. If not list of mutations
	# is supplied, then it plots the 10 mutations with the highest prevalence ratio.  The __init__ function took in a pre-established 
	# mutations dynamics dataframe that could be protein specific covariates, PANGO/WHO specific covariates, or a combination of
	# protein-pango, protein-who mutations; hence, to accurately label the visualization, the protein name or PANGO/WHO clade will need
	# to be passed in to the function.
	"""Parameters:
		mutations: Optional list of mutations to specificially plot
		protein: Optional protein name, but recomended if initialized dynamics dataframe is protein specific to  label plot
		PANGO: Optional Pango lineage name, but recomended if initialized dynamics dataframe is lineage specific to label plot
		WHO: Optional WHO name, but recommended if initialized dynamics dataframe is WHO name specific to label plot
	"""
	def plot_mutations(self, protein = "", PANGO = "", WHO = "", mutations = []):
		
		data = self.dynamics_df

		if WHO and PANGO:
			sys.exit("Input Error: WHO name and PANGO Lineage name not allowed to be iputted together.")

		if mutations:
			try:
				data = data[(data['Variant'].isin(mutations))]
			except:
				raise Exception("Invalid mutations input. Could not find all user supplied mutations.")

		if protein:
			data['Variant'] = data['Variant'].str.replace(protein+"_", "")

		top_mutations = data.sort_values(by = "Prevalence", ascending = False).drop_duplicates(subset = "Variant").head(15)
		top_mutations = list(top_mutations['Variant'])
		data = data[data['Variant'].isin(top_mutations)]
		plot_data = pd.DataFrame({'Variant': top_mutations})
		for date in self.dates:
			date_data = data[data['Date'] == date]
			year_month = date.to_period('M').strftime("%Y-%m")
			date_growth_data = pd.DataFrame({'Variant': list(date_data['Variant']), year_month: list(date_data['Prevalence'])})
			plot_data = plot_data.join(date_growth_data.set_index('Variant'), on='Variant')

		plot_data = plot_data.fillna(0)
		plot_data = plot_data.set_index('Variant')
		plot_data = plot_data[plot_data.columns[::-1]]
		x_val = np.asarray(plot_data.columns.tolist())
		indices = plot_data.index.tolist()
		for i in indices:
			mask = np.isfinite(plot_data.loc[i].values)
			plt.plot(x_val[mask], plot_data.loc[i].values[mask], label = i)
		plt.suptitle('A Time Course of Prevalence Ratios for '+PANGO+WHO+' '+protein+' Mutations from Sequences Isolated in '+self.region, fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Mutation Prevalence Ratio (Mutation Count/Isolates Count)')
		plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.10), ncol = 5, fancybox = True, shadow = True)
		plt.show()

