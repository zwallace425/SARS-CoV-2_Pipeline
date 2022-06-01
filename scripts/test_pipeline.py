import sys
from clean_mol_seq import CleanMolSeq as cms
from variant_dynamics import VariantDynamics as vd
from variant_analysis import VariantAnalysis as va


if __name__ == "__main__":

	reference = sys.argv[1]
	sequences = sys.argv[2]
	interval = sys.argv[3]
	world = sys.argv[4]

	if world == "T":
		world = True
	elif world == "F":
		world = False

	clean_sequences = "clean_sequences.fasta"

	kept_acc = cms.clean_mol_seqs(sequences, clean_sequences, 29400, 0.9999, True, 'SARS-CoV-2_Sequence_QC.json.gz')

	cov_prev, aa_prev, region_date_counts = vd.spatial_temporal_prevalence(reference, clean_sequences, qc_storage_file = 'SARS-CoV-2_Sequence_QC.json.gz')

	cov_prevalence_df = va.analyze_dynamics(cov_prev, region_date_counts, interval, world)
	aa_prevalence_df = va.analyze_dynamics(aa_prev, region_date_counts, interval, world)

	#aa_prevalence_df.to_csv("aa_prevalence_test.txt", sep = '\t', index = False)

	print(cov_prevalence_df)
	print('\n')
	print(aa_prevalence_df)



