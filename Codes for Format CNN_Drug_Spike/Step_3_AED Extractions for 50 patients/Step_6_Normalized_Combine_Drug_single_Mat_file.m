% Normalize dosage of all 50 patients based on the 95 percentile 
clc
clear all
close all

% This matlab file normalized the drug concentration with 95 percentile
% value no drug weights wonot be too small due to larger drug concentration
% values. 

% Input file location 
files = dir('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\Combined_Drug\*sid*');
% Save file location 
save_path = 'C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\Combined_Drug_Normalize\';
% Raw buffer to save 20 drug values among all patients 
RAW_Buff  = cell(1,20);
for files_ind = 1:length(files )
    file_ = files(files_ind);
    
    data = load(fullfile(file_.folder, file_.name));
    % store the drug data into a buffer 
    Drugs_=data.drugs_buffer;
    for i = 1:size(Drugs_,2)
        Drugs_Spec = Drugs_(:,i);
        Store_ = Drugs_Spec(Drugs_Spec~=0); %get all the non zero values of the drug
        % Store them in a buffer 
        RAW_Buff{1,i} =[ RAW_Buff{1,i};Store_];
    end
end
% Find 95CI for each Drug Dosage
%%
format long e
% store the 95 percentile value for each drug 
Drug_95_CI = zeros(1,20);
% Loop for each drug 
for j = 1:length(RAW_Buff)
    if isempty(RAW_Buff{1,j})
        % If drug data is not available then divide drug data with smallest
        % number 
        Drug_95_CI(1,j)=realmin;
    else
        % Drug data is available then store the 95 percentile of the drug 
        Drug_95_CI(1,j)=prctile(RAW_Buff{1,j},95);
    end
end
% Normalize the Drug value among all patients based on 95 percentile among all
% patients

for files_ind = 1:length(files )
    file_ = files(files_ind);
    % Load drug data
    data = load(fullfile(file_.folder, file_.name));
    Drugs_=data.drugs_buffer;
    Drug_info = data.Drug_info;
    % Divide the drug data with 95 percentile value. 
    Drugs_Normalize = Drugs_./Drug_95_CI;
    % Save file with normalized drug data. 
        save([save_path, file_.name],...
             'Drugs_Normalize', 'Drug_info')

end

