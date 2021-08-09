'''
This program will import the original NGS reads data of 100 different probes (50 in the CNSL region and 
50 in nonCNSL region) on different samples. Based on 4 known breakpoint positions (10-40, 20-40, 27-34, 32-38) 
of the CNSL gene region, this program characterize the samples of the potential genetic variation on this gene.
For each sample, it could be classified as either deletion or duplication of CNSL gene at one of the breakpoints,
or wildtype.
'''

import numpy as np
import pandas as pd
from utility import *

def main():

	cnsl_data = pd.read_csv('cnsl_data.csv.gz', index_col=[0])

	# normalization by sample
	# the average of reads across all probes within the nonCNSL region for each sample
	normalization_factor = np.mean(cnsl_data[list(cnsl_data)[51:]], axis = 1)

	cnsl_normalized = cnsl_data.copy()
	cnsl_normalized.iloc[:,1:] = cnsl_normalized.iloc[:,1:].div(normalization_factor, axis = 0)


	# normalization by probe
	normalization_factor_probe = np.mean(cnsl_normalized[list(cnsl_normalized)[1:]], axis = 0)
	plt.hist(normalization_factor_probe, bins = 50);

	cnsl_normalized_2 = cnsl_normalized.copy()

	cnsl_normalized_2.iloc[:,1:] = cnsl_normalized_2.iloc[:,1:].div(normalization_factor_probe, axis = 1)

	# recoganize "bad" probes
	# CNSL_probe_5, _23, _46
	cnsl_std = np.std(cnsl_normalized_2.iloc[:, 1:51], axis = 0) 

	bad_prob_ind = np.where(cnsl_std > 0.2)

	# drop the problematic probes
	cnsl_clean = cnsl_normalized_2.drop((cnsl_std[cnsl_std > 0.2]).index, axis = 1)


	annotation = []
	for i in range(cnsl_normalized_2.shape[0]):
	    annotation.append(variation_calling(cnsl_normalized_2, i))
	cnsl_clean["annotation"] = annotation

	annotation_by_ethnicity = ethnicity_variation.groupby(["ethnicity", "annotation"]).size().reset_index(name='counts')

	cnsl_clean.to_csv("cnsl_annotation.csv", index = False)
	annotation_by_ethnicity.to_csv("cnsl_annotation_by_ethnicity.csv", index = False)



# Call the main function to run the program
if __name__ == "__main__":
	main()
