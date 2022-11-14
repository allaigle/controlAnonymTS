# Analysis of an anonymization method

## analysis_anonym_meth1_part1.py

This python script allows to measure normalized multivariate dissimilarities between real patient timeseries and anonymized timeseries.

### Precondition : 
To have the following packages installed : 
sys, os, re, glob, time, random, statistics, csv, pandas, array, numpy, dtaidistance.
If necessary, please look at the **install.txt** file in the root folder. 

### Execution

Command line: **python3 analysis_anonym_meth1_part1.py path_pr path_pa path_analysis_anonym filename_csv**

The differents arguments:
1. **path_pr** corresponds to the path of the folder containing original multivariate timeseries.
2. **path_pa** corresponds to the path of the folder containing anonymized multivariate timeseries.
3. **path_analysis_anonym** corresponds to the folder in which the pdf graphical timeserie file will be saved. The script creates the directory if it doesn't exist.
4. **filename_csv** corresponds to the output filename of normalized dissimilarities. 

Example command line (from the root folder): "python3 tools/analysis_anonym_meth1_part1.py multivariate_original_dataset/ gener_simulated_data_meth1_prop-level_0.50_perturb-level_0.05/ analysis_anonym_meth1_prop-level_0.50_perturb-level_0.05/ distri_dissim_norm_meth1_prop-level_0.50_perturb-level_0.05.csv".
This command line will create statistical datasets for each physiological parameter.

## analysis_anonym_meth1_part2.py

This python script allows to generate the different csv files containing either average, standard deviation, median, minimun or maximum values from the distribution of normalized dissimilarities.

### Precondition : 
All packages cited in the analysis_anonym_meth1_part1.py script have been installed and are commun to both scripts (part1 and part2).

### Execution

Command line: **python3 analysis_anonym_meth1_part2.py path_pr path_pa path_analysis_anonym val_prop val_perturb**

The differents arguments:
1. **path_pr** corresponds to the path of the folder containing original multivariate timeseries.
2. **path_pa** corresponds to the path of the folder containing anonymized multivariate timeseries.
3. **path_analysis_anonym** corresponds to the path of the folder containing new csv files (e.g., dissim_norm files) and where the new analyses will be saved in.
4. **val_prop** corresponds to the proportion level that was applied to measure dissimilarities. 
5. **val_perturb** corresponds to the perturbation level that was applied to measure dissimilarities. 

Example command line (from the root folder): "python3 ./tools/analysis_anonym_meth1_part2.py ./multivariate_original_dataset/ ./gener_simulated_data_meth1_prop-level_0.50_perturb-level_0.05/ ./analysis_anonym_meth1_prop-level_0.50_perturb-level_0.05/ 0.50 0.05". 
This command line will create csv files containing statistical values (avg, std, med, min & max) for 0.50 proportion and 0.05 perturbation levels.

## analysis_anonym_meth1_part3.py

This python script allows to generate the different csv files containing statistical test results for one physiological parameter and a pair of parameters.

### Precondition : 

To have the following packages installed : 
itertools, scipy.
If necessary, please look at the **install.txt** file in the root folder. 

Note that some packages used in this script (analysis_anonym_meth1_part3) have already been installed for the two previous scripts (part1 and part2).

### Execution

Command line: **python3 analysis_anonym_meth1_part3.py path_analysis_anonym param_physio val_prop val_perturb**

The differents arguments:
1. **path_analysis_anonym** corresponds to the path of the folder containing new csv files (e.g., dissim_norm csv files) and where the new analyses will be saved in.
2. **param_physio** corresponds to the physiological parameter that want to be analyzed (FC, PAS, PAM, PAD)
3. **val_prop** corresponds to the proportion level that was applied to measure dissimilarities (e.g., 0.50).
4. **val_perturb** corresponds to the perturbation level that was applied to measure dissimilarities (e.g., 0.05).

Example command line (from the root folder): "python3 ./tools/analysis_anonym_meth1_part3.py ./analysis_anonym_meth1_prop-level_0.50_perturb-level_0.05/ FC 0.50 0.05". 
This command line will create csv files containing statistical values (avg, std, med, min & max) for 0.50 proportion and 0.05 perturbation levels and only for the FC physiological parameter.

## boxplot_meth1.R

This R script allows to generate and save the boxplot illustating normalized dissimilarity distributions between real and anonymized patients for all parameters (each pair of perturbation/proportion levels).

### Precondition : 
To have the following packages installed : 
readr, ggplot2.
For their installation, please look at the **install.txt** file in the root folder. 

### Execution

Command line: **Rscript ./tools/boxplot_meth1.R**

This command line will save the generated boxplot into the ./tools/ directory as a .png file.
