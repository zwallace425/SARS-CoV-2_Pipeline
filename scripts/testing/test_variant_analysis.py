import sys
import numpy as np
import pandas as pd
import variant_analysis as va


if __name__ == "__main__":

	prevalence_df = sys.argv[1]
	interval = sys.argv[2]
	world = sys.argv[3]

	if (sys.argv[3] == "T"):
		world = True
	elif (sys.argv[3] == "F"):
		world = False

	prevalence_df = pd.read_csv(prevalence_df, sep = '\t')

	dynamics = va.VariantAnalysis.analyze_dynamics(prevalence_df, interval, world)

	print(dynamics)
