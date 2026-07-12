%%% Medication Analysis %%%%

clc
clear all
close all

% This matlab file is used to calculate the drug concentration for the each
% patient for each drug based on the weight of the patients.

%%% Drug information %%%%%


Drug_info = ["levetiracetam","lacosamide","lorazepam","phenytoin",...
    "fosphenytoin","phenobarbital", "carbamazepine","valproate","divalproex","topiramate","clobazam",...
    "lamotrigine", "oxcarbazepine","diazepam","zonisamide","clonazepam","propofol","midazolam","ketamine","pentobarbital"];

Drug_info = Drug_info';

%%% Convert


% Load csv sheet for getting the patient age,gender and weight
[~,~, RAW_Age_Gender] = xlsread('MNR_SID_50_Header_New_Sid.xlsx');
RAW_Age_Gender = RAW_Age_Gender(2:end,:);
%%% Load the CNN label file for find Sid information
excelfiles = dir('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\CNN_Label\*.csv');

% excel sheet with patient weights in sean query if any patient weight is
% not available for a given data then match with the patient from the sean
% data
[RAW_Weight]  =readtable('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\CausalModeling_IIIC\data\bodyweights_second_query_Sean.xlsx','Sheet','Combined Data - Weight Only');
RAW_Weight = table2cell(RAW_Weight );

RAW_Weight_remove_header = RAW_Weight;

%%%% Below code is to Calculating weight for the patients
%%%% Calculating weight for the patients

for k =1:length(excelfiles )
    
    file = excelfiles(k);
    f = fullfile(file.folder,file.name);
    [~,~,RAW_HL]  =xlsread(f);
    %Data_Buff = zeros(size(RAW_HL,1),1);
    start_time_ = extractAfter(file.name,"_");
    start_time_CNN =extractBefore(start_time_,"_Bad");
    start_time_CNN = datetime(start_time_CNN,'InputFormat','yyyyMMdd_HHmmss');
    
    
    sid = cellstr(string(extractBefore(file.name,"_")));
    index  = find(strcmp(RAW_Age_Gender(1:end,3),  sid));
    
    %%%% Calculate weight for the each patient from the excel sheet and
    %%%% select nearest weight during the starting of EEG date.
    
    indx = find((cell2mat(RAW_Weight_remove_header(1:end,1)))== cell2mat(RAW_Age_Gender(index,1)));
    if isempty(indx)
        %%%% Weight did not find in the excel sheet
        Weight_patient(k) = NaN;
    else
        
        date_weight = RAW_Weight_remove_header(indx,2);
        
        
        for date_rec = 1: length(date_weight)
            diff_date(date_rec) =abs((date_weight{date_rec} - start_time_CNN));
        end
        
        %%%%%% finding the minimum distance between dates %%%%%%
        [Magi,Posi] = min(diff_date);
        
        
        %%% Coonvert weight to pounds to kg %%%%%
        Weight_patient(k) = str2num(RAW_Weight_remove_header{indx(Posi),3})*0.453592;
    end
    
    
    Gender(k) = RAW_Age_Gender( index ,4);
    Age(k) = RAW_Age_Gender{ index ,5};
    clear diff_date
    
end

%%% Convert weights which is less than 20 Kg to Nans
Weight_patient(Weight_patient < 20) = NaN;

%%% Append missing weights
nan_indx = find(isnan(Weight_patient));

%%% Non Nans
Non_nan_indx = find(~isnan(Weight_patient));
Gender_Not_Nan = Gender(1,Non_nan_indx);
Age_Not_Nan = Age(1,Non_nan_indx);
Weight_patient_Not_Nan = Weight_patient(1,Non_nan_indx);

% For a patients whose weights are missing calculate the mean weigth from
% the same gender and age cohort from the weights that are found.
% loop for each patient whose weight is missing
for i = 1: length(nan_indx)
    Gender_Weinan  =   Gender(nan_indx(1,i));
    % Find the indices of the gender based on the missing patient gender
    ind =  find(strcmp(  Gender_Not_Nan, Gender_Weinan ));
    % find the patients whose age are present in the given missing patient
    % gender
    Age_thres =Age_Not_Nan(ind);
    % find their weights too
    Weight_patient_gender  = Weight_patient_Not_Nan(ind);
    % find the means weights  for the patient whose weight is missing
    % Mean weight is calculated from the finding the cohort of patients
    % who are of same gender and age is wiht a rane og 10 years
    Weight_mean = mean(Weight_patient_gender (find(Age_thres>=(Age(nan_indx(1,i))-10) & Age_thres<=(Age(nan_indx(1,i))+10))));
    
    if isnan(Weight_mean)
        % If still canot find the weight then increase the range of age by 20 years
        Weight_mean = mean(Weight_patient_gender (find(Age_thres>=(Age(nan_indx(1,i))-20) & Age_thres<=(Age(nan_indx(1,i))+20))));
    end
    
    Weight_patient(nan_indx (i))=Weight_mean ;
    
    
end

%%%%%% Above code is to calculate weight %%%%%%%%%%

%%
% Once weight of the patient is calcualted then need to calculate the drug
% concentration. Based on the units we need to multiple weight of the
% patient that's the reason we calculate the weigths of the patient before .

for k =1:length(excelfiles )
    
    % Load the CNN data which is used for finding the start time of the EEG
    % and total length of the cEEG so drug data is extracted within the
    % cEEG availability
    file = excelfiles(k);
    f = fullfile(file.folder,file.name);
    %[~,~,RAW_HL]  =xlsread(f);
    RAW_HL  =readtable(f,'ReadVariableNames',false);
    RAW_HL = table2array(RAW_HL);
    
    % find the start time form the CN data
    start_time_ = extractAfter(file.name,"_");
    start_time_CNN =extractBefore(start_time_,"_Bad");
    start_time_CNN = datetime(start_time_CNN,'InputFormat','yyyyMMdd_HHmmss');
    
    
    
    
    %%%% Load Pre and post  epic medication data  calculated form step 1  %%%%%%
    folder_path = 'C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\Sid_Pre_Post_EPIC\';
    filename = sprintf('%s_MAR_Pre_Post_Epic.xlsx', extractBefore(file.name,"_"));
    f_epic = fullfile(folder_path,filename);
    [~,~,RAW]  =xlsread(f_epic);
    
    
    
    
    %%%%% Make directory to store medication %%%%%
    pname = 'C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\AED_2_Sec\' ;
    dname = sprintf('%s',extractBefore(file.name,"_"));
    mkdir(fullfile(pname,dname));
    
    % format time
    str = string(datestr( cell2mat(RAW(2:end,3)), 'HH:MM:SS' ));
    RAW(2:end,3) = cellstr(str);
    
    %%%% Remove header names %%%%
    RAW_Remove_Header = RAW(2:end,:);
    
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
    
    %%%%% finding drugs in the drug infor variable
    drug_find = unique_drugs(indx);
    
    % Loop for identify drugs
    
    % for ind_drug = 1: length(drug_find)
    
    
    for ind_drug = 1: length(drug_find)
        Date_repeat=cellstr(start_time_CNN + seconds((0:size(RAW_HL,1)-1)*2));
        % Seperate buffer NSAED and SAED information
        Data_Buff_SAED = zeros(size(RAW_HL,1),1)  ;
        Data_Buff_NSAED = zeros(size(RAW_HL,1),1)  ;
        
        % seperate buffer for saving the weighted NSAED and SAED data
        Data_Buff_Nor_Weight_SAED =  zeros(size(RAW_HL,1),1)  ;
        Data_Buff_Nor_Weight_NSAED =  zeros(size(RAW_HL,1),1)  ;
        %%%% finding the indecies of the each unique drug
        indx = find(strcmp(RAW_Remove_Header(:,7),drug_find(ind_drug)));
        
        %%% Columns needed %%%%%
        % 2 column of date
        % 3 column of time
        % 5 column dosage
        % 6 units
        % 7 drug name
        column_index = [2 3 5 6 7];
        Data_drug = RAW_Remove_Header(indx,column_index);
        % variable for storing the formated date and time
        Date_buff = [];
        Data_sort = [];
        for date_ind = 1: size(Data_drug,1)
            
            %%% Formating date and time for the particular drug and its all data %%%%
            Date_buff{date_ind,1} = datetime(Data_drug{date_ind,1},'InputFormat','MM/dd/yyyy') + timeofday( datetime(Data_drug{date_ind,2},'InputFormat','HH:mm:ss'));
            %%%% Variable to save datenum variable for sorting
            Data_sort(date_ind) = datenum(Date_buff{date_ind,1});
            
        end
        
        %%%% Sorting with respect to date and time
        [~,idx] = sort(Data_sort, 'ascend');
        
        %%% Sorting date variables
        Date_buff = Date_buff(idx,1);
        
        %%% Sorting data_drug matrix
        Data_drug = Data_drug(idx,:);
        %Data_Buff  = [];
        
        for seconds_ind = 1: length(Date_buff)
            % find start_time so can check if it with in the CNN start time
            %
            
            start_time  = (round(seconds(Date_buff{seconds_ind} -start_time_CNN)/2))+1;
            % check for Drug units is 'MG' , 'MG_PE','MG PE','MCG','MG/KG'
            
            if isequal(Data_drug{seconds_ind,4},'MG') ||  isequal(Data_drug{seconds_ind,4},'MG_PE') ||   isequal(Data_drug{seconds_ind,4},'MG PE') || isequal(Data_drug{seconds_ind,4},'MCG')  || isequal(Data_drug{seconds_ind,4},'MG/KG')
                
                % if the drug units are 'MG' , 'MG_PE','MG
                % PE','MCG','MG/KG' then drug is bolus anf infusion is
                % within a minute so end time is stat time + minute
                end_time = ( start_time + (60/2))-1;
                
                % Based on the units convert the drug into standard one.
                if isequal(Data_drug{seconds_ind,4},'MG') ||  isequal(Data_drug{seconds_ind,4},'MG_PE') ||  isequal(Data_drug{seconds_ind,4},'MG PE')
                    Drug_dosage = Data_drug{seconds_ind,3}*60;
                elseif isequal(Data_drug{seconds_ind,4},'MCG')
                    Drug_dosage = Data_drug{seconds_ind,3}*(60/1000);
                elseif  isequal(Data_drug{seconds_ind,4},'MG/KG')
                    Drug_dosage = Data_drug{seconds_ind,3}*60* Weight_patient(k) ;
                end
                
                % if the drug units are 'MG/HR , 'MG/KG/HR','MCG/KG/HR'
                % 'MCG/KG/MIN' then drug is continuous  so end time will
                % until the next srat time. 
                
            elseif isequal(Data_drug{seconds_ind,4},'MG/HR') ||  isequal(Data_drug{seconds_ind,4},'MG/KG/HR') ||   isequal(Data_drug{seconds_ind,4},'MCG/KG/HR') || isequal(Data_drug{seconds_ind,4},'MCG/KG/MIN')
                
                if seconds_ind==length(Date_buff)
                    end_time = length(Data_Buff_NSAED);
                else
                    %%%% Adding dosage until one step before end date so no need to add
                    %%%% +1 for end that is added for start
                    end_time =  (round(seconds(Date_buff{seconds_ind+1} - start_time_CNN)/2))+1;
                end
                
                % Based on the units convert the drug into standard one.
                if isequal(Data_drug{seconds_ind,4},'MG/HR')
                    Drug_dosage = Data_drug{seconds_ind,3};
                elseif   isequal(Data_drug{seconds_ind,4},'MG/KG/HR')
                    Drug_dosage = Data_drug{seconds_ind,3}* Weight_patient(k);
                elseif  isequal(Data_drug{seconds_ind,4},'MCG/KG/HR')
                    Drug_dosage = (Data_drug{seconds_ind,3}*Weight_patient(k))/1000;
                elseif  isequal(Data_drug{seconds_ind,4},'MCG/KG/MIN')
                    Drug_dosage = (Data_drug{seconds_ind,3}*60* Weight_patient(k))/1000;
                end
                
                
                
            end
            
            % check start time is max of 1 and start time it cannot be less
            % than 1 if so means data is not within cEEG range 
            
            % end time should be minimum of total length of cEE means it
            % cannot be over CEEG end time. 
      
            start_time = max(1,start_time);
            end_time = min(length(Data_Buff_NSAED),end_time)  ;
            if end_time > start_time
                
                % Save drug data based on the units 
                if isequal(Data_drug{seconds_ind,4},'MG') ||  isequal(Data_drug{seconds_ind,4},'MG_PE') ||   isequal(Data_drug{seconds_ind,4},'MG PE') || isequal(Data_drug{seconds_ind,4},'MCG')  || isequal(Data_drug{seconds_ind,4},'MG/KG')
                    
                    Data_Buff_NSAED (start_time:end_time) = Drug_dosage;
                    
                    % Normalized the drug data based on the weight of the
                    % patient 
                    Data_Buff_Nor_Weight_NSAED  (start_time:end_time) = Drug_dosage/Weight_patient(k);
                    
                elseif isequal(Data_drug{seconds_ind,4},'MG/HR') ||  isequal(Data_drug{seconds_ind,4},'MG/KG/HR') ||   isequal(Data_drug{seconds_ind,4},'MCG/KG/HR') || isequal(Data_drug{seconds_ind,4},'MCG/KG/MIN')
                    
                    Data_Buff_SAED (start_time:end_time) = Drug_dosage;
                    
                    Data_Buff_Nor_Weight_SAED  (start_time:end_time) = Drug_dosage/Weight_patient(k);
                end
                
            end
            
        end
        
        % Combine both NSAED and SAED information 
        Data_Buff = Data_Buff_NSAED + Data_Buff_SAED;
        Data_Buff_Nor_Weight = Data_Buff_Nor_Weight_NSAED + Data_Buff_Nor_Weight_SAED ;
        
        Data_Buff = num2cell(Data_Buff);
        Data_Buff_Nor_Weight = num2cell(Data_Buff_Nor_Weight);
        Data = [Date_repeat',Data_Buff,Data_Buff_Nor_Weight ];
        Data_Header = {'Date&Time','Drug Amount','Drug Amount Normalised Weight'};
        RAW_Buff_inc_header = [Data_Header ;Data];
        file_Name = extractBefore(file.name,"_");
        
        %%%% Saving exel file
        
        file_Name_cre = sprintf('%s_%s_2secWindow.mat',file_Name,drug_find{ind_drug});
        
        
        folder_name = ['C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\AED_2_Sec\',sprintf('%s',file_Name),'\'];
        folder=[folder_name   file_Name_cre];
        
        save(folder, 'RAW_Buff_inc_header');
    end
    
    
end