# Testing the object for loading gisaid metadata and getting variant counts

import pandas as pd
from datetime import datetime
from gisaid_metadata import GisaidMetadata

gisaid_data = pd.read_csv('metadata_tsv_2022_06_02/metadata.tsv', sep = '\t')

gisaid_data = gisaid_data[(gisaid_data['Collection date'] >= '2022-02-01')]

gisaid_data = gisaid_data.head(10000)

gm = GisaidMetadata(gisaid_data)

cov_prev, cov_region_dates = gm.covariate_counts()
print('\n')
print(cov_prev)
print('\n')
print(cov_region_dates)

aa_prev, aa_region_dates = gm.mutation_counts()
print('\n')
print(aa_prev)
print('\n')
print(aa_region_dates)

spike_prev, spike_region_dates = gm.spike_counts()
print('\n')
print(spike_prev)
print('\n')
print(spike_region_dates)
