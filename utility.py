'''
Utility functions supporting the del_dup_analysis program.
'''


import numpy as np
import pandas as pd

def variation_calling(data, n):
    signal_array = data.iloc[n, 1:51]
    result = check_breakpoint(signal_array, 10, 40)
    if result:
        return result + "_10"
    result = check_breakpoint(signal_array, 20, 40)
    if result:
        return result + "_20"    
    result = check_breakpoint(signal_array, 27, 34)
    if result:
        return result + "_27" 
    result = check_breakpoint(signal_array, 32, 38)
    if result:
        return result + "_32" 
    return "wt"
    
    
def check_breakpoint(signal_array, u, d):
    ind = np.setdiff1d(np.arange(u, d+1), bad_prob_ind)
    signal_mean = np.mean(signal_array[ind])
    if abs(signal_mean - 1.5) < 0.15:
        return "duplication"
    elif abs(signal_mean - 0.5) < 0.15:
        return "deletion"
    else:
        return
