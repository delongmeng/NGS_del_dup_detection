'''
Utility functions supporting the del_dup_analysis program.
'''


import numpy as np
import pandas as pd
import argparse


def get_input_args():
	'''
	This function parses the arguments from user input. 
	''' 

	# Create Parse using ArgumentParser
	parser = argparse.ArgumentParser()

	# Create command line arguments as mentioned above using add_argument() from ArguementParser method
	parser.add_argument('--input_path', type = str, default = 'cnsl_data.csv.gz',
	                help = 'Please define the input path (the default is: cnsl_data.csv.gz)')

	return parser.parse_args()   



def variation_calling(data, sample_id, bad_prob_ind):
	'''
	arguments: 
		data - dataframe, (we use the normalized data here; note: keep the bad probes in it)
		sample_id - int, row index of sample
		bad_prob_ind - a list of probe indexes that will be removed from analysis
	output:
		string, annotation of a certain sample
	'''
	signal_array = data.iloc[sample_id, 1:51]
	result = check_breakpoint(signal_array, 10, 40, bad_prob_ind)
	if result:
	    return result + "_10"
	result = check_breakpoint(signal_array, 20, 40, bad_prob_ind)
	if result:
	    return result + "_20"    
	result = check_breakpoint(signal_array, 27, 34, bad_prob_ind)
	if result:
	    return result + "_27" 
	result = check_breakpoint(signal_array, 32, 38, bad_prob_ind)
	if result:
	    return result + "_32" 
	return "wt"
    
    
def check_breakpoint(signal_array, u, d, bad_prob_ind):
	'''
	arguments: 
		signal_array - an array of the signal of all CNSL region probes
		u - int, the upstream bound of the breakpoint to check
		d - int, the downstream bound of the breakpoint to check
		bad_prob_ind - a list of probe indexes that will be removed from analysis
	output:
		a string of "duplication" or "deletion", or None if there's no hit
	'''
	ind = np.setdiff1d(np.arange(u, d+1), bad_prob_ind)
	signal_mean = np.mean(signal_array[ind])
	if abs(signal_mean - 1.5) < 0.15:
	    return "duplication"
	elif abs(signal_mean - 0.5) < 0.15:
	    return "deletion"
	else:
	    return