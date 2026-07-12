CNN_Label : Formatted label data for the 50 patients,  once CNN applied to cEEG it label cEEG for every 2 second window. This 2 second data is concatenated
            so we will have a continuous data here NaNs are labelled for missing data or discontinue data. 
Combined_Drug_Normalize: Formatted Drug data for the 50 patients, Drug concentration is calculated for every 2 seconds based on the each patient weight 
			and normalized with maximum value of concentration
Spike Data : Formatted data for 50 patients. formatted output of spike detector into 2 seconds window and concatenated and NaNs are included when data is missing. 

Sim_Results_Store: Is the folder which has the results that are shown in the paper. 

Code_To_Generate_Figures_Paper: Consist of all codes that are used to generate figures for the paper. 

Step_1_Calculate_parameters: Is the matlab file which will calculate lognormal and drug weights for each patients. Since we are looking for 48 patients 
			      data it will generate excel sheet with parameters for 48 patients. 

Curve_fit_Python : Python code for curve fit function which include lognormal and drug concentration weights . This function is automically called inside 
		   Step_1_Calculate_parameters matlab file. 

Step_2_Calculate_Sample_Effect_Size: Matlab file for calculating the sample size and effect size using the parameters from the 48 patients. Input for the file is 
                                     Treatment group population select for complete cessation , outcomes and delay in drug response and treatment delay. 

Drug_Outcomes: Is a matlab file which is a function for calculating the outcomes for Step_2_Calculate_Sample_Effect_Size matlab file.