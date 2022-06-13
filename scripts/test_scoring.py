# Script for testing variant_scoring.py

import sys
import pandas as pd
import pickle
from variant_analysis import VariantAnalysis as va
from variant_scoring import VariantScoring as vs

with open('data/gisaid_cov_prevalence_df.pkl', 'rb') as handle:
	cov_prevalence_df = pickle.load(handle)

with open('data/gisaid_aa_prevalence_df.pkl', 'rb') as handle:
	aa_prevalence_df = pickle.load(handle)

with open('data/gisaid_spike_prevalence_df.pkl', 'rb') as handle:
	spike_prevalence_df = pickle.load(handle)

with open('data/gisaid_cov_region_date_counts_df.pkl', 'rb') as handle:
	cov_region_dates_df = pickle.load(handle)

with open('data/gisaid_aa_region_date_counts_df.pkl', 'rb') as handle:
	aa_region_dates_df = pickle.load(handle)

with open('data/gisaid_spike_region_date_counts_df.pkl', 'rb') as handle:
	spike_region_dates_df = pickle.load(handle)

spike_sfoc = pd.read_csv("data/Spike_SFoCs.txt", sep = '\t')
non_spike_sfoc = pd.read_csv("data/Non-Spike_SFoCs.txt", sep = '\t')

cov_prevalence_df = va.analyze_dynamics(cov_prevalence_df, cov_region_dates_df, interval = 'M')
aa_prevalence_df = va.analyze_dynamics(aa_prevalence_df, aa_region_dates_df, interval = 'M')
spike_prevalence_df = va.analyze_dynamics(spike_prevalence_df, spike_region_dates_df, interval = 'M')

cov = vs(cov_prevalence_df)
aa = vs(aa_prevalence_df)
spike = vs(spike_prevalence_df)

#print("Sequence Prevalence Score: Full")
#print(cov.sequence_prevalence_score().head(25))
#print('\n')

#print("Sequence Prevalence Score: Spike")
#print(spike.sequence_prevalence_score().head(25))
#print('\n')

#print("Mutation Prevalence Score: Spike")
#print(aa.mutation_prevalence_score(Spike = True).head(25))
#print('\n')

#print("Mutation Prevalence Score: Non-Spike")
#print(aa.mutation_prevalence_score(nonSpike = True).head(25))
#print('\n')

print("Composite Score: Full")
print(cov.composite_score(spike_sfoc = spike_sfoc, non_spike_sfoc = non_spike_sfoc).head(25))
print('\n')

print("Composite Score: Full, Spike")
print(cov.composite_score(spike_sfoc = spike_sfoc).head(25))
print('\n')

print("Composite Score: Full, Non-Spike")
print(cov.composite_score(non_spike_sfoc = non_spike_sfoc).head(25))
print('\n')

print("Composite Score: Spike")
print(spike.composite_score(spike_sfoc = spike_sfoc).head(25))
print('\n')


