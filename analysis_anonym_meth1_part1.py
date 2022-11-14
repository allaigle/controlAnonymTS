#!/usr/bin/env python3

###############################################################################
# Programm allowing to produce "distri_dissim_norm_meth1[_par_val].csv" files #
#         Author: Alice Laigle - M2BB - Nantes University (2021/2022)         #
###############################################################################

# Import libraries
import sys, os, re, glob, time
import random, statistics
import csv, pandas as pd
import array, numpy as np
from dtaidistance import dtw


""" Usage of arguments in bash command line """
path_pr = sys.argv[1] # path of the original folder
path_pa = sys.argv[2] # path of the anonymized folder
path_analysis_anonym = sys.argv[3] # path of the analysis_anonym_meth1_<param1>_<value1>_<param2>_<value2> folder 
filename_csv = sys.argv[4] # output : distri_dissim_norm_meth1_<param1>_<value1>_<param2>_<value2>.csv

########################################################################################################################################################

def get_random_pr_series(path_pr):
    """ Return a list of 10 random files from the real "series" files
    Argument:
    path_real = path of the original folder (real patients, multivariate) """
    
    try:
        list_pr = glob.glob(os.path.join(path_pr, "*.txt")) # Check files into the original folder 
        list_pr_series = [] # Create an empty list that will only contains "series" files
      
        for files_pr in list_pr: 
            if files_pr.endswith("_series.txt"):   # Keep only files that are finishing by "_series.txt"
                list_pr_series.append(files_pr)  # If True, add the file to the list 
            
        list_pr_series.sort(key=lambda f: int(re.sub('\D', '', f))) # Sort the series' list by ascending order with lambda key & regular expression library   

        random_pr_series = random.sample(list_pr_series, 10) # Pick 10 random files from the list (rand.seed() is placed in the "main" part)

    except FileNotFoundError:
        sys.stderr.write(f"[FileNotFoundError] Impossible to open the folder : {path_pr} \n")
        exit(1)
    
    return random_pr_series

########################################################################################################################################################

def get_list_pa_series(path_pa):
    """Return the list of 'series' files for anonymized patient 
    Argument:
    path_pa = path of anonymized folder """
    
    try:
        list_pa_multi = glob.glob(os.path.join(path_pa, "*.txt")) # Check files into the anonymized folder 
        list_pa_series = []
     
        for files_pa in list_pa_multi: 
            if files_pa.endswith("_series.txt"):
                list_pa_series.append(files_pa)
            
    except FileNotFoundError:
        sys.stderr.write(f"[FileNotFoundError] Impossible to open the folder : {path_pa} \n")
        exit(1)
        
    return list_pa_series

#######################################################################################################################################################

def DTWu_m_minimum(path_pr, path_pa, path_analysis_anonym):
    """ Read series files, measure univariate DTW distances between each anonymized patients and 10 real patients 
    (randomly chosen), then measure multivariate DTW disantces (mean of FC, PAS, PAM & PAD), and write in a  
    temporary file all minimum DTWm
    Arguments:
    path_pr = path of the original folder def natural_sort(l): 
    path_pa = path of anonymized folder 
    path_analysis_anonym = path of the analysis_anonym_meth1 folder 
    """
    
    list_pa_series = get_list_pa_series(path_pa) # Get anonymized patient series
    list_pa_series.sort(key=lambda f: int(re.sub('\D', '', f))) # Sort the series' list by ascending order      
    lists_dist_param, list_arrays, list_DTWm, DTWm_minimum = [], [], [], [] # Create lists for to store data of nested loops 
    
    # Create a file in the analysis_anonym_meth1 folder to temporary memorize minimum DTWm 
    temp_DTWm_minimum = pd.DataFrame() # First, create a dataframe that will contains minimum DTWm
    temp_DTWm_minimum.to_csv(f'{path_analysis_anonym}/temp_DTWm_minimum.txt', sep = ',', header=False, index=False) # Create the file (updated in next loop)
   
    # For the 1000 anonymized patients :
    for pa in list_pa_series :
        
        list_dist_10_fc, list_10_dist_pas, list_10_dist_pam, list_10_dist_pad = [], [], [], [] # Need local variables here (not outside of the loop) 
                                                                                                # in order to re-initialize them for each anonymized patient
        # Read anonymized time series and transform pd.DataFrame into np.ndarray
        pa_series_fc = (pd.read_csv(pa, usecols = [1])).to_numpy() # 1 : FC column (where [0] is for Time)
        pa_series_pas = (pd.read_csv(pa, usecols = [2])).to_numpy() # 2 : PAS column 
        pa_series_pam = (pd.read_csv(pa, usecols = [3])).to_numpy() # 3 : PAM column 
        pa_series_pad = (pd.read_csv(pa, usecols = [4])).to_numpy() # 4 : PAD column 
        
        random_pr_10_series = get_random_pr_series(path_pr) # Get 10 random pr series for each pa series
        
        for i in random_pr_10_series :
            # Read real time series
            pr_series_fc = (pd.read_csv(i, usecols = [1])).to_numpy() # 1 : FC column 
            pr_series_pas = (pd.read_csv(i, usecols = [2])).to_numpy() # 2 : PAS column  
            pr_series_pam = (pd.read_csv(i, usecols = [3])).to_numpy() # 3 : PAM column 
            pr_series_pad = (pd.read_csv(i, usecols = [4])).to_numpy() # 4 : PAD column 

            # Measure univariate DTW distances between pa and pr
            distance_fc = dtw.distance(array.array('d', pa_series_fc), array.array('d', pr_series_fc), use_c = True) # pa versus pr_1 (pr_1+1 for each turn)
            list_dist_10_fc.append(distance_fc) # Add the distance to the list 
            distance_pas = dtw.distance(array.array('d', pa_series_pas), array.array('d', pr_series_pas), use_c = True) 
            list_10_dist_pas.append(distance_pas)
            distance_pam = dtw.distance(array.array('d', pa_series_pam), array.array('d', pr_series_pam), use_c = True) 
            list_10_dist_pam.append(distance_pam)
            distance_pad = dtw.distance(array.array('d', pa_series_pad), array.array('d', pr_series_pad), use_c = True) 
            list_10_dist_pad.append(distance_pad) 
        
        # Measure multivariate DTW 
        lists_dist_param = [list_dist_10_fc, list_10_dist_pas, list_10_dist_pam, list_10_dist_pad] # Make a list of the four lists
        list_arrays = [np.array(x) for x in lists_dist_param] # For each sublist, make np.arrays (e.g., [[2, 3], [4, 6]] --> to then compare 2 with 4 & 3 with 6)
        list_DTWm = [np.mean(k) for k in zip(*list_arrays)] # Calculate the mean of each parameter distances (e.g., mean_of_dist_1)    
        
        # Memorize the minimum DTWm on the 10 real patients 
        list_DTWm = sorted(list_DTWm, key = lambda x:float(x)) # Sort the list of strings in an ascendent order
        DTWm_minimum.append(list_DTWm[0]) # Keep the first element of each sorted list_DTWm 
        if ((len(DTWm_minimum)) == (len(list_pa_series))) :  # Once all minimum DTWm have been calculated, write them in the temp file 
            with open(f'{path_analysis_anonym}/temp_DTWm_minimum.txt', 'a') as f:    # 'a' is for appending, quite optional here 
                writer = csv.writer(f)      
                writer.writerow(DTWm_minimum)
    
    return

########################################################################################################################################################  

def normalize_DTWm_min(path_analysis_anonym):
    
    """ Calculate the mean E and the standard error S from the temporary DTWm_minimum file 
    and return each normalized dissimilarities (1000) into a csv file
    Argument : 
    path_analysis_anonym = path of the analysis_anonym_meth1 folder  
    """

    try :
        # Find and read the temporary DTWm_minimum file
        csv_file = glob.glob(os.path.join(path_analysis_anonym, "temp_DTWm_minimum.txt"))
        for f in csv_file:
            DTWm_minimum = (pd.read_csv(f)).columns.tolist() # Read the file and convert the type 'pandas.core.frame.DataFrame' into a list
            
    except IsADirectoryError:
        sys.stderr.write(f"[IsADirectoryError] Is a directory: '{path_analysis_anonym}' \n")
        exit(1)
        
    except FileNotFoundError :
        sys.stderr.write(f"[FileNotFoundError] No such file or directory: 'temp_DTWm_minimum.txt', {path_analysis_anonym} \n")
        exit(1)  
        
    try : 
        DTWm_minimum = [float(items) for items in DTWm_minimum] # Convert the DTWm_minimum list of strings into a list of floats
        mean_E = statistics.mean(DTWm_minimum) # Calculate the mean of all minimum DTWm
        standev_S = statistics.stdev(DTWm_minimum) # Calculate the standard deviation of all minimum DTWm
        DTWm_minimum_array = np.array(DTWm_minimum) # Convert the list of float into a numpy array
        dissim_norm = (DTWm_minimum_array - mean_E)/standev_S # Calculate normalized dissimilarities
        
    except ValueError :
        sys.stderr.write(f"[StatisticError] There is a problem in at least one statistics operation.\n")
        exit(1)
    
    try :    
        # Save normalized dissimilarities as a csv file (txt works too) with header 'dissim_norm' and delimiter ','
        np.savetxt(f'{path_analysis_anonym}{filename_csv}', dissim_norm, delimiter = ',', header = "dissim_norm", fmt='%f', comments = '')
 
    except IOError :
        sys.stderr.write(f"[IOError] No such file or directory : {path_analysis_anonym} \n")
        exit(1)

########################################################################################################################################################  

def main():
    
    random.seed(time.time()) # Generate seeds depending on time (as random as possible)
    
    print('BEGIN : Calculate the 1000 univariate, multivariate DTW distances and keep minimum DTWm distances into a temporary file. \n')
    print('Note : This script takes around 1m20.')
    DTWu_m_minimum(path_pr, path_pa, path_analysis_anonym) # Run the DTWu_m_minimum function

    print('BEGIN : Calculate the distribution of the 1000 minimum DTWm. \n') 
    normalize_DTWm_min(path_analysis_anonym) # Rn the normalize_DTWm_min function
    
    print('REMOVE : temporary file containing the 1000 minimum DTWm. \n')
    os.remove(f'{path_analysis_anonym}/temp_DTWm_minimum.txt') # Suppress the temporary file from the analysis_anonym_meth1 folder
    
    print('END OF : analysis_anonym_meth1_part1.py. \n') 
    
if __name__ == '__main__' : 
    main()