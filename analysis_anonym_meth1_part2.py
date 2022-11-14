#!/usr/bin/env python3

########################################################################
# Programm producing statistical values for each parameter of datasets #
#      Author: Alice Laigle - M2BB - Nantes University (2021/2022)     #
########################################################################

# Import libraries
import sys, os, re, glob
import pandas as pd
import numpy as np
import statistics


""" Usage of arguments in bash command line """
path_pr = sys.argv[1] # path of the original folder
path_pa = sys.argv[2] # path of the anonymized folder
path_analysis_anonym = sys.argv[3] # path of the analysis_anonymized folder
val_prop = sys.argv[4] # value of the proportional level tested
val_perturb = sys.argv[5] # value of the perturbation level tested

#####################################################################################################################################################

def get_list_series(path_pr, path_pa):
    """Return lists of 'series' files for real and anonymized patients 
    Argument:
    path_pr = path of original folder
    path_pa = path of anonymized folder """
    
    try:
        # Check files into the original and anonymized folders 
        list_pr_multi = glob.glob(os.path.join(path_pr, "*.txt")) # Original folder 
        list_pa_multi = glob.glob(os.path.join(path_pa, "*.txt")) # Anonymized folder 
        list_pr_series, list_pa_series = [], [] # Creation of empty lists to add series' files inside

        for files_pa in list_pr_multi:     # For real patients 
            if files_pa.endswith("_series.txt"):
                list_pr_series.append(files_pa)
     
        for files_pa in list_pa_multi:      # For anonymized ones
            if files_pa.endswith("_series.txt"):
                list_pa_series.append(files_pa)
        
        # Sort series' lists by ascending order   
        list_pa_series.sort(key = lambda f: int(re.sub('\D', '', f)))  
        list_pr_series.sort(key = lambda f: int(re.sub('\D', '', f))) 
            
    except FileNotFoundError:
        sys.stderr.write(f"[FileNotFoundError] Impossible to open at least one of the folders : {path_pr} and/or {path_pa} \n")
        exit(1)
        
    return list_pr_series, list_pa_series


#######################################################################################################################################################

def calculate_values(path_pr, path_pa, path_analysis_anonym, val_prop, val_perturb): 
    
    # Read anonymized time series and transform pd.DataFrame into np.ndarray
    list_pr_series, list_pa_series = get_list_series(path_pr, path_pa) # Get series' lists* 
    
    avg_pa_fc, avg_pa_pas, avg_pa_pam, avg_pa_pad, std_pa_fc, std_pa_pas, std_pa_pam, std_pa_pad = ([] for i in range(8))
    med_pa_fc, med_pa_pas, med_pa_pam, med_pa_pad = ([] for i in range(4))
    min_pa_fc, min_pa_pas, min_pa_pam, min_pa_pad, max_pa_fc, max_pa_pas, max_pa_pam, max_pa_pad = ([] for i in range(8)) 
    
    for pa in list_pa_series :
        
        pa_series = pd.read_csv(pa)
        
        pa_series_fc = pa_series['FC'].values.tolist() # Convert the FC into a list of values         
        pa_series_pas, pa_series_pam, pa_series_pad = pa_series['PAS'].values.tolist(), pa_series['PAM'].values.tolist(), pa_series['PAD'].values.tolist()  

        # Calculte means and stdev for each pa file, add it to a np.array and convert it to a pd.df
        avg_pa_fc, std_pa_fc = pd.DataFrame(np.append(avg_pa_fc, (statistics.mean(pa_series_fc)))), pd.DataFrame(np.append(std_pa_fc, (statistics.stdev(pa_series_fc))))  
        avg_pa_pas, std_pa_pas = pd.DataFrame(np.append(avg_pa_pas, (statistics.mean(pa_series_pas)))), pd.DataFrame(np.append(std_pa_pas, (statistics.stdev(pa_series_pas))))  
        avg_pa_pam, std_pa_pam = pd.DataFrame(np.append(avg_pa_pam, (statistics.mean(pa_series_pam)))), pd.DataFrame(np.append(std_pa_pam, (statistics.stdev(pa_series_pam))))  
        avg_pa_pad, std_pa_pad = pd.DataFrame(np.append(avg_pa_pad, (statistics.mean(pa_series_pad)))), pd.DataFrame(np.append(std_pa_pad, (statistics.stdev(pa_series_pad))))  

        # Calculte median for each pa file, add it to a np.array and convert it to a pd.df
        med_pa_fc = pd.DataFrame(np.append(med_pa_fc, (statistics.median(pa_series_fc))))
        med_pa_pas = pd.DataFrame(np.append(med_pa_pas, (statistics.median(pa_series_pas))))
        med_pa_pam = pd.DataFrame(np.append(med_pa_pam, (statistics.median(pa_series_pam))))
        med_pa_pad = pd.DataFrame(np.append(med_pa_pad, (statistics.median(pa_series_pad))))
       
        # Sort series_param by ascending order
        sort_pa_fc, sort_pa_pas = sorted(pa_series_fc, key = float), sorted(pa_series_pas, key = float)
        sort_pa_pam, sort_pa_pad  = sorted(pa_series_pam, key = float), sorted(pa_series_pad, key = float) 
        
        # MIN
        min_pa_fc = pd.DataFrame(np.append(min_pa_fc, (sort_pa_fc[0]))) # keep the first element (= minimum value)
        min_pa_pas = pd.DataFrame(np.append(min_pa_pas, (sort_pa_pas[0]))) 
        min_pa_pam = pd.DataFrame(np.append(min_pa_pam, (sort_pa_pam[0]))) 
        min_pa_pad = pd.DataFrame(np.append(min_pa_pad, (sort_pa_pad[0]))) 
        # MAX
        max_pa_fc = pd.DataFrame(np.append(max_pa_fc, (sort_pa_fc[-1]))) # keep the last element (= maximim value)
        max_pa_pas = pd.DataFrame(np.append(max_pa_pas, (sort_pa_pas[-1]))) 
        max_pa_pam = pd.DataFrame(np.append(max_pa_pam, (sort_pa_pam[-1]))) 
        max_pa_pad = pd.DataFrame(np.append(max_pa_pad, (sort_pa_pad[-1])))

    # For real patients : 
    avg_pr_fc, avg_pr_pas, avg_pr_pam, avg_pr_pad, std_pr_fc, std_pr_pas, std_pr_pam, std_pr_pad = ([] for i in range(8))
    med_pr_fc, med_pr_pas, med_pr_pam, med_pr_pad = ([] for i in range(4))
    min_pr_fc, min_pr_pas, min_pr_pam, min_pr_pad, max_pr_fc, max_pr_pas, max_pr_pam, max_pr_pad = ([] for i in range(8)) 
    
    for pr in list_pr_series :
        
        pr_series = pd.read_csv(pr)
        
        pr_series_fc = pr_series['FC'].values.tolist() # Convert the FC into a list of values         
        pr_series_pas, pr_series_pam, pr_series_pad = pr_series['PAS'].values.tolist(), pr_series['PAM'].values.tolist(), pr_series['PAD'].values.tolist()  

        # Calculte means and stdev for each pr file, add it to a np.array and convert it to a pd.df
        avg_pr_fc, std_pr_fc = pd.DataFrame(np.append(avg_pr_fc, (statistics.mean(pr_series_fc)))), pd.DataFrame(np.append(std_pr_fc, (statistics.stdev(pr_series_fc))))  
        avg_pr_pas, std_pr_pas = pd.DataFrame(np.append(avg_pr_pas, (statistics.mean(pr_series_pas)))), pd.DataFrame(np.append(std_pr_pas, (statistics.stdev(pr_series_pas))))  
        avg_pr_pam, std_pr_pam = pd.DataFrame(np.append(avg_pr_pam, (statistics.mean(pr_series_pam)))), pd.DataFrame(np.append(std_pr_pam, (statistics.stdev(pr_series_pam))))  
        avg_pr_pad, std_pr_pad = pd.DataFrame(np.append(avg_pr_pad, (statistics.mean(pr_series_pad)))), pd.DataFrame(np.append(std_pr_pad, (statistics.stdev(pr_series_pad))))  

        # Calculte median for each pa file, add it to a np.array and convert it to a pd.df
        med_pr_fc = pd.DataFrame(np.append(med_pr_fc, (statistics.median(pr_series_fc))))
        med_pr_pas = pd.DataFrame(np.append(med_pr_pas, (statistics.median(pr_series_pas))))
        med_pr_pam = pd.DataFrame(np.append(med_pr_pam, (statistics.median(pr_series_pam))))
        med_pr_pad = pd.DataFrame(np.append(med_pr_pad, (statistics.median(pr_series_pad))))
       
        # Sort series_param by ascending order
        sort_pr_fc, sort_pr_pas = sorted(pr_series_fc, key = float), sorted(pr_series_pas, key = float)
        sort_pr_pam, sort_pr_pad  = sorted(pr_series_pam, key = float), sorted(pr_series_pad, key = float) 
        
        # Calculate the minimum for each original file from sorted lists, add it to a np.array and convert it to a pd.df
        min_pr_fc = pd.DataFrame(np.append(min_pr_fc, (sort_pr_fc[0]))) # [0] : keep the first element (= minimum value)
        min_pr_pas = pd.DataFrame(np.append(min_pr_pas, (sort_pr_pas[0]))) 
        min_pr_pam = pd.DataFrame(np.append(min_pr_pam, (sort_pr_pam[0]))) 
        min_pr_pad = pd.DataFrame(np.append(min_pr_pad, (sort_pr_pad[0]))) 
        # Calculate the maximum for each original file from sorted lists, add it to a np.array and convert it to a pd.df
        max_pr_fc = pd.DataFrame(np.append(max_pr_fc, (sort_pr_fc[-1]))) # [-1] : keeps the last element (= maximim value)
        max_pr_pas = pd.DataFrame(np.append(max_pr_pas, (sort_pr_pas[-1]))) 
        max_pr_pam = pd.DataFrame(np.append(max_pr_pam, (sort_pr_pam[-1]))) 
        max_pr_pad = pd.DataFrame(np.append(max_pr_pad, (sort_pr_pad[-1])))

    # Save files : AVG # distri_dissim_norm_meth1_prop-level_0.25_perturb-level_0.05
    avg_fc = pd.merge(avg_pa_fc, avg_pr_fc, left_index=True, right_index=True, suffixes=('avg_anonym', 'avg_real')) # Merge both anonimyzed and real mean dataframes for FC 
    avg_fc.to_csv(f'{path_analysis_anonym}/avg_values_meth1_FC_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['avg_anonym', 'avg_real'], sep = ',', index = False) # Save it as a csv file  
    avg_pas = pd.merge(avg_pa_pas, avg_pr_pas, left_index=True, right_index=True, suffixes=('avg_anonym', 'avg_real')) 
    avg_pas.to_csv(f'{path_analysis_anonym}/avg_values_meth1_PAS_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['avg_anonym', 'avg_real'], sep = ',', index = False) 
    avg_pam = pd.merge(avg_pa_pam, avg_pr_pam, left_index=True, right_index=True, suffixes=('avg_anonym', 'avg_real')) 
    avg_pam.to_csv(f'{path_analysis_anonym}/avg_values_meth1_PAM_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['avg_anonym', 'avg_real'], sep = ',', index = False) 
    avg_pad = pd.merge(avg_pa_pad, avg_pr_pad, left_index=True, right_index=True, suffixes=('avg_anonym', 'avg_real')) 
    avg_pad.to_csv(f'{path_analysis_anonym}/avg_values_meth1_PAD_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['avg_anonym', 'avg_real'], sep = ',', index = False) 
    # STD
    std_fc = pd.merge(std_pa_fc, std_pr_fc, left_index=True, right_index=True, suffixes=('std_anonym', 'std_real'))
    std_fc.to_csv(f'{path_analysis_anonym}/std_values_meth1_FC_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['std_anonym', 'std_real'], sep = ',', index = False) 
    std_pas = pd.merge(std_pa_pas, std_pr_pas, left_index=True, right_index=True, suffixes=('std_anonym', 'std_real')) 
    std_pas.to_csv(f'{path_analysis_anonym}/std_values_meth1_PAS_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['std_anonym', 'std_real'], sep = ',', index = False) 
    std_pam = pd.merge(std_pa_pam, std_pr_pam, left_index=True, right_index=True, suffixes=('std_anonym', 'std_real')) 
    std_pam.to_csv(f'{path_analysis_anonym}/std_values_meth1_PAM_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['std_anonym', 'std_real'], sep = ',', index = False) 
    std_pad = pd.merge(std_pa_pad, std_pr_pad, left_index=True, right_index=True, suffixes=('std_anonym', 'std_real')) 
    std_pad.to_csv(f'{path_analysis_anonym}/std_values_meth1_PAD_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['std_anonym', 'std_real'], sep = ',', index = False) 
    # MED
    med_fc = pd.merge(med_pa_fc, med_pr_fc, left_index=True, right_index=True, suffixes=('med_anonym', 'med_real')) 
    med_fc.to_csv(f'{path_analysis_anonym}/med_values_meth1_FC_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['med_anonym', 'med_real'], sep = ',', index = False)   
    med_pas = pd.merge(med_pa_pas, med_pr_pas, left_index=True, right_index=True, suffixes=('med_anonym', 'med_real')) 
    med_pas.to_csv(f'{path_analysis_anonym}/med_values_meth1_PAS_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['med_anonym', 'med_real'], sep = ',', index = False) 
    med_pam = pd.merge(med_pa_pam, med_pr_pam, left_index=True, right_index=True, suffixes=('med_anonym', 'med_real')) 
    med_pam.to_csv(f'{path_analysis_anonym}/med_values_meth1_PAM_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['med_anonym', 'med_real'], sep = ',', index = False) 
    med_pad = pd.merge(med_pa_pad, med_pr_pad, left_index=True, right_index=True, suffixes=('med_anonym', 'med_real')) 
    med_pad.to_csv(f'{path_analysis_anonym}/med_values_meth1_PAD_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['med_anonym', 'med_real'], sep = ',', index = False) 
    # MIN
    min_fc = pd.merge(min_pa_fc, min_pr_fc, left_index=True, right_index=True, suffixes=('min_anonym', 'min_real')) 
    min_fc.to_csv(f'{path_analysis_anonym}/min_values_meth1_FC_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['min_anonym', 'min_real'], sep = ',', index = False) 
    min_pas = pd.merge(min_pa_pas, min_pr_pas, left_index=True, right_index=True, suffixes=('min_anonym', 'min_real')) 
    min_pas.to_csv(f'{path_analysis_anonym}/min_values_meth1_PAS_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['min_anonym', 'min_real'], sep = ',', index = False) 
    min_pam = pd.merge(min_pa_pam, min_pr_pam, left_index=True, right_index=True, suffixes=('min_anonym', 'min_real')) 
    min_pam.to_csv(f'{path_analysis_anonym}/min_values_meth1_PAM_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['min_anonym', 'min_real'], sep = ',', index = False) 
    min_pad = pd.merge(min_pa_pad, min_pr_pad, left_index=True, right_index=True, suffixes=('min_anonym', 'min_real')) 
    min_pad.to_csv(f'{path_analysis_anonym}/min_values_meth1_PAD_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['min_anonym', 'min_real'], sep = ',', index = False) 
    # MAX
    max_fc = pd.merge(max_pa_fc, max_pr_fc, left_index=True, right_index=True, suffixes=('max_anonym', 'max_real')) # Merge both anonimyzed and real mean dataframes for FC 
    max_fc.to_csv(f'{path_analysis_anonym}/max_values_meth1_FC_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['max_anonym', 'max_real'], sep = ',', index = False) # Save it as a csv file  
    max_pas = pd.merge(max_pa_pas, max_pr_pas, left_index=True, right_index=True, suffixes=('max_anonym', 'max_real')) 
    max_pas.to_csv(f'{path_analysis_anonym}/max_values_meth1_PAS_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['max_anonym', 'max_real'], sep = ',', index = False) 
    max_pam = pd.merge(max_pa_pam, max_pr_pam, left_index=True, right_index=True, suffixes=('max_anonym', 'max_real')) 
    max_pam.to_csv(f'{path_analysis_anonym}/max_values_meth1_PAM_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['max_anonym', 'max_real'], sep = ',', index = False) 
    max_pad = pd.merge(max_pa_pad, max_pr_pad, left_index=True, right_index=True, suffixes=('max_anonym', 'max_real')) 
    max_pad.to_csv(f'{path_analysis_anonym}/max_values_meth1_PAD_prop-level_{val_prop}_perturb-level_{val_perturb}.csv', header = ['max_anonym', 'max_real'], sep = ',', index = False) 
       
#######################################################################################################################################################  

def main():
    
    print('\nBEGIN : Calculate means, standard deviations, medians, minimum and maximum values for anonymized and real patients.')
    print('Note : This script takes around 15 seconds. \n')
    
    calculate_values(path_pr, path_pa, path_analysis_anonym, val_prop, val_perturb)
    
    print('END OF : analysis_anonym_meth1_part2.py. \n') 
    
if __name__ == '__main__' : 
    main()