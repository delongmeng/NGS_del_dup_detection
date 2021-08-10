'''
This program will import the original NGS reads data of 100 different probes (50 in the CNSL region and 
50 in nonCNSL region) on different samples. Based on 4 known breakpoint positions (10-40, 20-40, 27-34, 32-38) 
of the CNSL gene region, this program characterize the samples of the potential genetic variation on this gene.
For each sample, it could be classified as either deletion or duplication of CNSL gene at one of the breakpoints,
or wildtype.
'''

import numpy as np
import pandas as pd
import argparse
from utility import *

def main():

	# take arguments from the command (the input path and output path)
	in_arg = get_input_args()

	# load data
	cnsl_data = pd.read_csv(in_arg.input_path, index_col=[0])

	# sample normalization 
	# using the average of reads across all probes within the nonCNSL region for each sample
	sample_normalization_factor = np.mean(cnsl_data[list(cnsl_data)[51:]], axis = 1)
	cnsl_s_norm = cnsl_data.copy()
	cnsl_s_norm.iloc[:,1:] = cnsl_s_norm.iloc[:,1:].div(sample_normalization_factor, axis = 0)

	# probe normalization
	# using the average of reads across all samples within the nonCNSL region for each probe
	probe_normalization_factor = np.mean(cnsl_s_norm[list(cnsl_s_norm)[1:]], axis = 0)
	cnsl_s_p_norm = cnsl_s_norm.copy()
	cnsl_s_p_norm.iloc[:,1:] = cnsl_s_p_norm.iloc[:,1:].div(probe_normalization_factor, axis = 1)

	# recoganize "bad" probes
	cnsl_std = np.std(cnsl_s_p_norm.iloc[:, 1:51], axis = 0) 
	bad_prob_ind = np.where(cnsl_std > 0.2)

	# drop the problematic probes
	cnsl_clean = cnsl_s_p_norm.drop((cnsl_std[cnsl_std > 0.2]).index, axis = 1)

	# preform annotation
	annotation = []
	for i in range(cnsl_s_p_norm.shape[0]):
	    annotation.append(variation_calling(cnsl_s_p_norm, i, bad_prob_ind))
	cnsl_clean["annotation"] = annotation

	# summarize the annotations according to ethnicity groups
	ethnicity_variation = cnsl_clean[["ethnicity", "annotation"]]
	annotation_by_ethnicity = ethnicity_variation.groupby(["ethnicity", "annotation"]).size().reset_index(name='counts')
	ethnicity_counts = cnsl_data.ethnicity.value_counts().reset_index(name='total_counts').rename(columns = {"index": "ethnicity"})
	annotation_by_ethnicity = pd.merge(annotation_by_ethnicity, ethnicity_counts, how = "left", on = "ethnicity")
	annotation_by_ethnicity["frequency"] = annotation_by_ethnicity.counts/annotation_by_ethnicity.total_counts
	breakpoint_by_ethnicity = annotation_by_ethnicity.replace(["deletion", "duplication"], "breakpoint", regex = True)
	breakpoint_by_ethnicity = breakpoint_by_ethnicity.groupby(["ethnicity","annotation"]).agg({'counts':'sum','total_counts':'mean', 'frequency':'sum'}).reset_index()

	# save output files
	cnsl_clean.to_csv("cnsl_annotation.csv", index = False)
	annotation_by_ethnicity.to_csv("cnsl_annotation_by_ethnicity.csv", index = False)
	breakpoint_by_ethnicity.to_csv("cnsl_breakpoint_by_ethnicity.csv", index = False)


# Call the main function to run the program
if __name__ == "__main__":
	main()
