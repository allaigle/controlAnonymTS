#!/usr/bin/env python3

################################################################
# Programm producing statistical tests from statistical values #
#   for one physiological parameter and a pair of parameters    #
#  Author: Alice Laigle - M2BB - Nantes University (2021/2022) #
################################################################

# Import libraries
import sys, os, glob
import random, time
import pandas as pd
import itertools
from scipy.stats import ks_2samp, mannwhitneyu

""" Usage of arguments in bash command line """ 
path_analysis_anonym = sys.argv[1] # path of the analysis_anonymized folder
param_physio = sys.argv[2] # physiological parameter that want to be analyzed (FC, PAS, PAM, PAD)
val_prop = sys.argv[3] # value of the proportion level (e.g., 0.50)
val_perturb = sys.argv[4] # value of the perturbation level (e.g., 0.25)


# Generate the 'tests_meth1_<param_physio>_avg_<par>_<valeur>.csv x4 (one per param_phy) x4 (one per pair of prop-/perturb-levels)

#####################################################################################################################################################

def read_stat_csv(path_analysis_anonym, param_physio, val_prop, val_perturb):
    """ Return five dataframes of read lists (one for each statistical value)
    Arguments:
    path_analysis_anonym = path of the analysis_anonym_meth1_<par1>_<val1>_<par2>_<val2> folder that will be tested 
    param_physio = one of the four physiological parameters that will be tested (FC, PAS, PAM or PAD)
    val_prop = value of the proportion level that will be tested
    val_perturb = alue of the perturbation level that will be tested """
    
    try:
        # Make the list of all csv files from the analysis_anonym_meth1_<par1>_<val1>_<par2>_<val2> folder
        stat_list = glob.glob(os.path.join(path_analysis_anonym,"*.csv"))
        list_avg, list_std, list_med, list_min, list_max  = ([] for i in range(5))

        # Create lists of the relevant files
        for files in stat_list : 
            if files.endswith("avg_values_meth1_{}_prop-level_{}_perturb-level_{}.csv".format(param_physio, val_prop, val_perturb)): 
                list_avg.append(files)
            elif files.endswith("std_values_meth1_{}_prop-level_{}_perturb-level_{}.csv".format(param_physio, val_prop, val_perturb)):
                list_std.append(files)
            elif files.endswith("med_values_meth1_{}_prop-level_{}_perturb-level_{}.csv".format(param_physio, val_prop, val_perturb)):
                list_med.append(files)
            elif files.endswith("min_values_meth1_{}_prop-level_{}_perturb-level_{}.csv".format(param_physio, val_prop, val_perturb)):
                list_min.append(files)
            elif files.endswith("max_values_meth1_{}_prop-level_{}_perturb-level_{}.csv".format(param_physio, val_prop, val_perturb)):
                list_max.append(files)
            # Else would be either the dissim_norm file or files with other physiolocal parameters.

        # Read each list
        for a, s, me, mi, ma in itertools.product(list_avg, list_std, list_med, list_min, list_max): 
            df_read_avg = pd.read_csv(a, usecols=['avg_anonym', 'avg_real']) 
            df_read_std = pd.read_csv(s, usecols=['std_anonym', 'std_real'])
            df_read_med = pd.read_csv(me, usecols=['med_anonym', 'med_real'])
            df_read_min = pd.read_csv(mi, usecols=['min_anonym', 'min_real'])
            df_read_max = pd.read_csv(ma, usecols=['max_anonym', 'max_real'])     

    except FileNotFoundError:
        sys.stderr.write(f"[FileNotFoundError] No such file or directory. \n")
        exit(1)
        
    return df_read_avg, df_read_std, df_read_med, df_read_min, df_read_max 


#######################################################################################################################################################

def tests_KS_WMW_up():
    """ Perform Kolmogorow-Smirnov and Mann-Whitney U tests for each statistical value (avg, std, med, min & max)
    from the given physiological and valued parameters (passed as arguments of the command line)
    Note : our anonymization method is an unpaired method
    """
    try :
        # Read files 
        df_read_avg, df_read_std, df_read_med, df_read_min, df_read_max = read_stat_csv(path_analysis_anonym, param_physio, val_prop, val_perturb) 
    
        ## AVG 
        arr_avg_anonym, arr_avg_real = df_read_avg['avg_anonym'].to_numpy(), df_read_avg['avg_real'].to_numpy() # Transform the dataframe into 2 arrays to make tests
        avg_statKS, avg_pvalKS = ks_2samp(arr_avg_anonym, arr_avg_real) # Perform KS test on both samples and keep stat and pvalue into variables
        avg_statWMW, avg_pvalWMW = mannwhitneyu(arr_avg_anonym, arr_avg_real) # Perform WMW test on both samples and keep stat and pvalue into variables
        # Place result into dataframes
        avg_KS, avg_WMW = pd.DataFrame({'statKS': [avg_statKS], 'pvalKS': [avg_pvalKS]}), pd.DataFrame({'statWMW_up':[avg_statWMW], 'pvalWMW_up' : [avg_pvalWMW]})
        avg_result = pd.concat([avg_KS, avg_WMW], axis = 1, join = 'inner') # Concatenate both dataframes
        # Save these results as a csv file with headers : statKS, pvalKS, statWMW_up, pvalWMW_up and delimiters as commas
        avg_result.to_csv(f'{path_analysis_anonym}/test_meth1_{param_physio}_avg_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = True, sep = ',', index = False)
  
        ## STD 
        arr_std_anonym, arr_std_real = df_read_std['std_anonym'].to_numpy(), df_read_std['std_real'].to_numpy()
        std_statKS, std_pvalKS = ks_2samp(arr_std_anonym, arr_std_real) 
        std_statWMW, std_pvalWMW = mannwhitneyu(arr_std_anonym, arr_std_real) 
        std_KS, std_WMW = pd.DataFrame({'statKS': [std_statKS], 'pvalKS': [std_pvalKS]}), pd.DataFrame({'statWMW_up':[std_statWMW], 'pvalWMW_up' : [std_pvalWMW]})
        std_result = pd.concat([std_KS, std_WMW], axis = 1, join = 'inner') 
        std_result.to_csv(f'{path_analysis_anonym}/test_meth1_{param_physio}_std_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = True, sep = ',', index = False)
        
        ## MED 
        arr_med_anonym, arr_med_real = df_read_med['med_anonym'].to_numpy(), df_read_med['med_real'].to_numpy()
        med_statKS, med_pvalKS = ks_2samp(arr_med_anonym, arr_med_real) 
        med_statWMW, med_pvalWMW = mannwhitneyu(arr_med_anonym, arr_med_real) 
        med_KS, med_WMW = pd.DataFrame({'statKS': [med_statKS], 'pvalKS': [med_pvalKS]}), pd.DataFrame({'statWMW_up':[med_statWMW], 'pvalWMW_up' : [med_pvalWMW]})
        med_result = pd.concat([med_KS, med_WMW], axis = 1, join = 'inner') 
        med_result.to_csv(f'{path_analysis_anonym}/test_meth1_{param_physio}_med_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = True, sep = ',', index = False)

        ## MIN 
        arr_min_anonym, arr_min_real = df_read_min['min_anonym'].to_numpy(), df_read_min['min_real'].to_numpy()
        min_statKS, min_pvalKS = ks_2samp(arr_min_anonym, arr_min_real) 
        min_statWMW, min_pvalWMW = mannwhitneyu(arr_min_anonym, arr_min_real) 
        min_KS, min_WMW = pd.DataFrame({'statKS': [min_statKS], 'pvalKS': [min_pvalKS]}), pd.DataFrame({'statWMW_up':[min_statWMW], 'pvalWMW_up' : [min_pvalWMW]})
        min_result = pd.concat([min_KS, min_WMW], axis = 1, join = 'inner') 
        min_result.to_csv(f'{path_analysis_anonym}/test_meth1_{param_physio}_min_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = True, sep = ',', index = False)

        ## MAX 
        arr_max_anonym, arr_max_real = df_read_max['max_anonym'].to_numpy(), df_read_max['max_real'].to_numpy()
        max_statKS, max_pvalKS = ks_2samp(arr_max_anonym, arr_max_real) 
        max_statWMW, max_pvalWMW = mannwhitneyu(arr_max_anonym, arr_max_real) 
        max_KS, max_WMW = pd.DataFrame({'statKS': [max_statKS], 'pvalKS': [max_pvalKS]}), pd.DataFrame({'statWMW_up':[max_statWMW], 'pvalWMW_up' : [max_pvalWMW]})
        max_result = pd.concat([max_KS, max_WMW], axis = 1, join = 'inner') 
        max_result.to_csv(f'{path_analysis_anonym}/test_meth1_{param_physio}_max_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = True, sep = ',', index = False)

    except TypeError :
        sys.stderr.write(f"[TypeError] Error when trying to merge dataframes or wrong type of passed arguments. \n")
        exit(1)
    
    except ValueError:
        sys.stderr.write(f"[ValueError] NaN values might be present in some dataframes. \n")
        exit(1)

#######################################################################################################################################################  

def main():
    
    random.seed(time.time()) # Generate seeds depending on time (as random as possible)
    
    print('\nBEGIN : Calculate statKS, pvalKS, statWMW_up, pvalWMW_up and save them into a csv file for each physiological parameter and statistic values (avg, std, med, min, max).')
    tests_KS_WMW_up()
    
    print('END OF : analysis_anonym_meth1_part3.py. \n') 
    
if __name__ == '__main__' : 
    main()