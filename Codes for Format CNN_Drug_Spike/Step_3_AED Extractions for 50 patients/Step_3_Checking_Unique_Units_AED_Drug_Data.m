%%% Medication Analysis %%%%

clc
clear all
close all

% This matlab file is used to extract all unique units of the drugs so if
% there is any new units can find how to convert them into drug
% concentration.

%%% Drug information %%%%%
Drug_info = ["levetiracetam","lacosamide","lorazepam","phenytoin",...
    "fosphenytoin","phenobarbital", "carbamazepine","valproate","divalproex","topiramate","clobazam",...
    "lamotrigine", "oxcarbazepine","diazepam","zonisamide","clonazepam","propofol","midazolam","ketamine","pentobarbital"];


Drug_info = Drug_info';

% Load all drug information extracted from the Pre post epic for the each
% patient
excelfiles = dir('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\Sid_Pre_Post_EPIC\*.xlsx');

% buffer to save units
unique_units = {};

% load each patient pre post peic excel file.
for k =1:length(excelfiles )
    file = excelfiles(k);
    f = fullfile(file.folder,file.name);
    [~,~,RAW]  =xlsread(f);
    file_Name = extractBefore(file.name,"_");
    
    
    
    % format time.
    str = string(datestr( cell2mat(RAW(2:end,3)), 'HH:MM:SS' ));
    RAW(2:end,3) = cellstr(str);
    
    %%%% Remove header names %%%%
    RAW_Remove_Header = RAW(2:end,:);
    
    % Combining the drug which are similar
    %%% Combining phenytoin = fosphenytoin  &&&& valproate = divalproex
    indx_phy = find(strcmp(RAW_Remove_Header(:,7),"phenytoin"));
    
    for i_phy = 1:length(indx_phy)
        RAW_Remove_Header{indx_phy(i_phy),7} =  'fosphenytoin';
    end
    
    %%% Combining   valproate = divalproex
    indx_dival = find(strcmp(RAW_Remove_Header(:,7),"divalproex"));
    
    for i_dival = 1:length(indx_dival)
        RAW_Remove_Header{indx_dival(i_dival),7} =  'valproate';
    end
    
    
    
    %%%% Finding unique drugs
    unique_drugs = unique(cellfun(@num2str,RAW_Remove_Header(:,7),'UniformOutput',false));
    %%%% Find if unique drug is present in the drug info
    indx = ismember(unique_drugs,Drug_info);
    
    %%%%% finding drugs in the drug info variable
    drug_find = unique_drugs(indx);
    
    % Loop for identify drugs
    
    for ind_drug = 1: length(drug_find)
        
        %%%% finding the indecies of the each unique drug
        indx = find(strcmp(RAW_Remove_Header(:,7),drug_find(ind_drug)));
        
        %%% Columns needed %%%%%
        
        % column 6 is drug units
        % column 7 is drug name
        
        unique_units_buff = unique(cellfun(@num2str,RAW_Remove_Header(indx,6),'UniformOutput',false));
        unique_units = [unique_units;unique_units_buff];
    end
    unique_units_final = unique(cellfun(@num2str,unique_units,'UniformOutput',false));
end

%%% Store unique drug units in the excel file
file_Name = sprintf('Unique_Units_Sida.xlsx');
folder=['C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\'  file_Name];

Data_Header(1,1) ="Units";
Unique_Units = [Data_Header ;unique_units_final];
xlswrite(folder, Unique_Units);




