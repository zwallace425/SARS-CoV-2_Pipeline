# Object for scoring the variants with the BV-BRC designed scoring hueristics, including sequence/mutation prevalence
# score, functional impact score, and composite score.  This only accepts a dataframe in the format returned
# from analyze_dynamics in the VariantAnalysis object, which is ...

"""	Variant 		Region		Date 		Variant Count 		Isolates Count 		Prevalence 		Growth 		Jerk
	X				X			YYYY-MM-DD	X 					X 					X 				X 			X
"""

import sys
import re
import pandas as pd
from datetime import datetime

class VariantScoring(object):

	# Intialize the variant dynamics dataframe to only contain the last three timestamps of data
	def __init__(self, variants_df, interval = 3):
		
		dates = variants_df.sort_values(by = 'Date', ascending = False).drop_duplicates(subset = ['Date']).head(interval)
		dates = list(dates['Date'])
		dynamics_df = variants_df[(variants_df['Date'].isin(dates))]
		
		self.variants = dynamics_df


	# This function computes the composite score for a covariate, combining both the Sequence Prevalence Score and
	# the Functional Impact Score.
	"""Parameters:
		spike_sfoc: Spike protein Sequence Featurs of Concern dataframe curated by the BV-BRC team
		non_spike_sfocs: Non-Spike proteins Sequences Features of Concern (for now, nsp3, nsp5, nsp12, ORF3a, and N)
			curated by the BV-BRC team
		covariates: Optional list of covariates to specifically compute the composite score for
	"""
	def composite_score(self, spike_sfoc = [], non_spike_sfoc = [], covariates = []):

		print("Computing the Composite Score ...")

		if (len(spike_sfoc) == 0 and len(non_spike_sfoc) == 0):
			sys.exit("Input Error: Must input Spike SFoC dataframe, Non-Spike SFoC dataframe, or both.")
		
		prevalence_score = VariantScoring.sequence_prevalence_score(self, covariates = covariates)

		if covariates:
			impact_score = VariantScoring.functional_impact_score(covariates = covariates, spike_sfoc = spike_sfoc, non_spike_sfoc = non_spike_sfoc)
		else:
			impact_score = VariantScoring.functional_impact_score(covariates = list(prevalence_score['Variant']), spike_sfoc = spike_sfoc, non_spike_sfoc = non_spike_sfoc)
		
		composite_score = pd.merge(prevalence_score, impact_score, on = "Variant", how = "outer").fillna(0)
		composite_score['Composite Score'] = composite_score['Sequence Prevalence Score'] + composite_score['Functional Impact Score']
		composite_ranking = composite_score.sort_values(by = ['Composite Score'], ascending = [False]).drop_duplicates().reset_index(drop=True).dropna()

		print("Done computing the Composite Score")
		print('\n')
		
		return(composite_ranking.reset_index(drop=True))


	# This function computes the sequence prevalence scores of variants and return the variants ranked by score
	"""Parameters:
		covariates: Optional list of covariates to specifically compute the sequences prevalence score for
	"""
	def sequence_prevalence_score(self, covariates = []):
		
		data = self.variants

		if (covariates):
			data = data[(data['Variant'].isin(covariates))]

		significant_variants = data[(data['Variant Count'] > 10) & (data['Date'] == max(data['Date']))]
		#significant_variants = data
		significant_variants = significant_variants[['Variant', 'Region']]
		significant_variants_df = significant_variants.merge(data, on = ['Variant', 'Region']).drop_duplicates()
		significant_variants_df = significant_variants_df[(significant_variants_df['Prevalence'] > 0.05) | (significant_variants_df['Growth'] > 5)]

		scores = significant_variants_df.groupby('Variant').size().reset_index(
			name = 'Sequence Prevalence Score').sort_values(by = 'Sequence Prevalence Score', ascending = False).dropna()


		return scores.reset_index(drop=True)

	# This function computes the mutation prevalence score of AA mutations and returns the variants ranked by score.
	# It takes in an optional boolean parameter to nonSpike and if it's true it will only compute mutation prevalence
	# scoring for nonSpike proteins.
	"""Parameters:
		Spike: Boolean True or false to determine if only Spike mutations should be scored
		nonSpike: Boolean True or false to determine if only nonSpike mutations should be scored
	"""
	def mutation_prevalence_score(self, Spike = False, nonSpike = False):

		print("Computing Mutation Prevalence Score ...")
		
		if (Spike and nonSpike):
			sys.exit("Invalid Input: Only Spike or nonSpike may be set to boolean True, got True and True")

		data = self.variants

		if Spike:
			data = data[(data['Variant'].str.contains('Spike'))]
		
		if nonSpike:
			data = data[(data['Variant'].str.contains('Spike') == False)]

		significant_variants = data[(data['Variant Count'] > 10) & (data['Date'] == max(data['Date']))]
		significant_variants = significant_variants[['Variant', 'Region']]
		significant_variants_df = significant_variants.merge(data, on = ['Variant', 'Region']).drop_duplicates()
		significant_variants_df = significant_variants_df[(significant_variants_df['Prevalence'] > 0.05) | (significant_variants_df['Growth'] > 5)]

		scores = significant_variants_df.groupby('Variant').size().reset_index(
			name = 'Mutation Prevalence Score').sort_values(by = 'Mutation Prevalence Score', ascending = False).dropna()

		print("Done computing Mutation Prevalence Score.")
		print('\n')

		return scores.reset_index(drop=True)


	# This function compute the emerging lineage score for PANGO Lineages and returns the lineages ranked by the score.
	# This function does not recieve any parameters. NOTE: There must be a PANGO Lineage column with the intialized 
	# dataframe mapping each covariate to the PANGO Lineage or else the function will abort.
	def emerging_lineage_score(self):

		print("Computing Emerging Lineage Score ...")

		data = self.variants

		if 'PANGO Lineage' not in data.columns:
			sys.exit("Data Error: PANGO Lineage names not in the data --- cannot compute Emerging Lineage Score")

		significant_variants = data[(data['Variant Count'] > 10) & (data['Date'] == max(data['Date']))]
		significant_variants = significant_variants[['Variant', 'Region']]
		significant_variants_df = significant_variants.merge(data, on = ['Variant', 'Region']).drop_duplicates()
		significant_variants_df = significant_variants_df[(significant_variants_df['Growth'] > 15)]

		scores = significant_variants_df.groupby('PANGO Lineage').size().reset_index(
			name = 'Emerging Lineage Score').sort_values(by = 'Emerging Lineage Score', ascending = False).dropna()

		print("Done computing Emerging Lineage Score.")
		print('\n')

		return scores.reset_index(drop=True)


	# This function computes the functional impact score for variants within the last three dates of data.  If a spike
	# Sequence Features of Concern dataframe is inputted, it will compute Spike functional impact.  If a non-Spike
	# Sequence Features of Concern dataframe is inputted, it will compute non-Spike functional impact.  If both 
	# are inputted, it will compute a combined Spike and non-Spike functional impact scoring.  If neither, the 
	# function aborts.
	"""Parameters:
		covariates: A list of covariates to compute functional impact score. Note, there should be no duplicates
		spike_sfoc: Spike protein Sequence Featurs of Concern dataframe curated by the BV-BRC team
		non_spike_sfocs: Non-Spike proteins Sequences Features of Concern (for now, nsp3, nsp5, nsp12, ORF3a, and N)
			curated by the BV-BRC team
	"""
	@staticmethod
	def functional_impact_score(covariates, spike_sfoc = [], non_spike_sfoc = []):
		
		if (len(spike_sfoc) == 0 and len(non_spike_sfoc) == 0):
			sys.exit("Input Error: Must input Spike SFoC dataframe, Non-Spike SFoC dataframe, or both.")

		spike_impact = []
		non_spike_impact = []

		if (len(spike_sfoc) > 0):
			scores = []
			for i in range(len(covariates)):
				scores.append(VariantScoring.spike_functional_score(covariates[i], spike_sfoc))
			
			spike_impact = pd.DataFrame({'Variant': covariates, 'Spike Functional Score': scores})
			spike_impact = spike_impact.sort_values(by = 'Spike Functional Score', ascending = False)

		if (len(non_spike_sfoc) > 0):
			scores = []
			for i in range(len(covariates)):
				scores.append(VariantScoring.non_spike_functional_score(covariates[i], non_spike_sfoc))
			
			non_spike_impact = pd.DataFrame({'Variant': covariates, 'Non-Spike Functional Score': scores})
			non_spike_impact = non_spike_impact.sort_values(by = 'Non-Spike Functional Score', ascending = False)

		if (len(spike_impact) > 0 and len(non_spike_impact) > 0):
			functional_impact = spike_impact.merge(non_spike_impact, on = "Variant")
			functional_impact['Functional Impact Score'] = functional_impact['Spike Functional Score'] + functional_impact['Non-Spike Functional Score']
		elif (len(spike_impact) > 0 and len(non_spike_impact) == 0):
			functional_impact = spike_impact.rename({'Spike Functional Score': 'Functional Impact Score'}, axis = 1)
		elif (len(spike_impact) == 0 and len(non_spike_impact) > 0):
			functional_impact = non_spike_impact.rename({'Non-Spike Functional Score': 'Functional Impact Score'}, axis = 1)

		return functional_impact.reset_index(drop=True)


	# This is a helper function meant for computing the Spike functional impact score.  It will compute the Spike
	# functional impact score for a single covariate.  The supplied sfoc_df must be the Spike Sequence Features of
	# Concern dataframe
	"""Parameters:
		sfoc_df: Spike Sequence Features of Concern dataframe curated by the BV-BRC team
	"""
	@staticmethod
	def spike_functional_score(covariate, sfoc_df):
		mutations = covariate.split(",")
		sfoc_score = 0
		prev_pos = 1
		prev_mut = "N"
		for j in range(len(mutations)):
			mutation = mutations[j].split("_")
			if "Spike" not in mutation[0]:
				continue
			aa_pos = int(re.sub("[^0-9]", "", mutation[1]))
			aa_muts = re.sub(r'[0-9]+', '', mutation[1])
			aa_muts = aa_muts.replace("del", "-")
			mut = aa_muts[len(aa_muts) - 1]
			sfocs = sfoc_df[(sfoc_df['Start'] <= aa_pos) & (sfoc_df['End'] >= aa_pos)]
			sfocs.reset_index(inplace = True, drop = True)
			if (not ((aa_pos == (prev_pos + 1)) and (mut == "-") and (prev_mut	== "-"))):
				for k in range(len(sfocs)):
					if (aa_pos == 614):
						sfoc_score += 1
					if (pd.notna(sfocs.at[k,"mAb escape"])):
						if ("class 1" in sfocs.loc[k].at["mAb escape"]):
							sfoc_score += 1
						if ("class 2" in sfocs.loc[k].at["mAb escape"]):
							sfoc_score += 1
						if ("class 3" in sfocs.loc[k].at["mAb escape"]):
							sfoc_score += 1
						if ("class 4" in sfocs.loc[k].at["mAb escape"]):
							sfoc_score += 1
					if (pd.notna(sfocs.at[k,"serum Ab escape"])):
						if ("convalescent serum" in sfocs.loc[k].at["serum Ab escape"]):
							sfoc_score += 1
						if ("Moderna vaccine serum" in sfocs.loc[k].at["serum Ab escape"]):
							sfoc_score += 1
					if (pd.notna(sfocs.at[k,"Increased ACE2 binding"])):
						sfoc_score += 1
					if (pd.notna(sfocs.at[k,"Region of interest"])):
						sfoc_score += 1
			prev_pos = aa_pos
			prev_mut = mut

		return sfoc_score


	# This is a helper function meant for computing the non-Spike functional impact score.  It will compute the non-Spike
	# functional impact score for a single covariate.  The supplied sfoc_df must be the non-Spike Sequence Features of
	# Concern dataframe
	"""Parameters:
		sfoc_df: Non-Spike Sequence Features of Concern dataframe curated by the BV-BRC team
	"""
	@staticmethod
	def non_spike_functional_score(covariate, sfoc_df):
		mutations = covariate.split(",")
		sfoc_score = 0
		for j in range(len(mutations)):
			mutation = mutations[j].split("_")
			if "Spike" in mutation[0]:
				continue
			prot = mutation[0]
			aa_pos = int(re.sub("[^0-9]", "", mutation[1]))
			aa_muts = re.sub(r'[0-9]+', '', mutation[1])
			mut = aa_muts[len(aa_muts) - 1]
			sfocs = sfoc_df[(sfoc_df['Protein'] == prot) & (sfoc_df['Start'] <= aa_pos) & (sfoc_df['End'] >= aa_pos)]
			sfocs.reset_index(inplace = True, drop = True)
			sfoc_score += len(sfocs)
			
		return sfoc_score
		









