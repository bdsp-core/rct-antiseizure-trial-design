clc
clear all
close all
% AED_2_Sec folder is having the individual drug data for each patient this
% matlab code is used to combine all individual drug data into an one buffer.
% Combined drug data has 20 columns and rows equal to size of CNN labels.
% If a drug is available for the patient then the particular column of the
% drug is over lap with the actual drug concentration of the patient or else the column will be zeros.

Drug_info = ["levetiracetam","lacosamide","lorazepam","phenytoin",...
    "fosphenytoin","phenobarbital", "carbamazepine","valproate","divalproex","topiramate","clobazam",...
    "lamotrigine", "oxcarbazepine","diazepam","zonisamide","clonazepam","propofol","midazolam","ketamine","pentobarbital"];
% Input file with all drug infomration for each client
files = dir('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\AED_2_Sec\*sid*');
% save path
save_path = 'C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\Combined_Drug\'

% load drug concentration data for each patient
for files_ind = 1:length(files )
    
    file_ = files(files_ind);
    sid_file = dir(fullfile(file_.folder, file_.name,'*.mat'));
    % find length of the mat file so can generate a mat of size 20X length
    data = load(fullfile(sid_file(1).folder, sid_file(1).name));
    % rows in drug buffer is equal to rows in sids drug data
    drugs_buffer = zeros(size(data.RAW_Buff_inc_header,1)-1,20);
    
    % Loop for each drug for a patient loop
    for sid_file_ind = 1:length(sid_file)
        
        % look for drug present in the Sid
        if ismember(extractBetween(sid_file(sid_file_ind).name,"_","_"),Drug_info)
            % find location of drug in the drug info for assigning in
            % the buffer
            [Sid_logic,drug_loc] = ismember(extractBetween(sid_file(sid_file_ind).name,"_","_"),Drug_info);
            % load drug data
            data_drug = load(fullfile(sid_file(sid_file_ind).folder, sid_file(sid_file_ind).name));
            data_drug_ = data_drug.RAW_Buff_inc_header;
            % taking only the normalized data
            data_drug_  = cell2mat(data_drug_(2:end,3));
            % Assign the weight normalized data to the buffer based on
            % drug location
            drugs_buffer(:,drug_loc)=data_drug_;
            
        end
        
    end
    % save the data
    save([save_path, file_.name],...
        'drugs_buffer', 'Drug_info')
    
    
end
