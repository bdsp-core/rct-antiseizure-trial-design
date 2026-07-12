%%% Data Analysis MAR Pre Epic Data %%%%

clc
clear all
close all



% Load xlsx files Pre epic medication sheet saved in local D folder
excelfiles = dir('D:\eMAR_Processed\*.xlsx');
% Load MNR SID excel file 
%%% first coloumn is mnr and second column is Sid
%Data_MNR_Sid = {};
[~,~,Data_MNR_Sid] = xlsread('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Info 50 clients in SAH SAGE DataBase.xlsx');

Data_MRN_Header = Data_MNR_Sid(1,:);
% Remove header 
Data_MNR_Sid = Data_MNR_Sid(2:end,:);


% Cell array Buffer to save each patient medication data for combining both pre and
% post epic and post and pre epic data is concatenated 
RAW_Buff  = cell(size(Data_MNR_Sid,1),1);

%%% Loading post epic medication data 
[~,~,RAW_Post] = xlsread('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\CausalModeling_IIIC\data\all-meds-from-list-20180118.csv');

%%%% Remove header 

RAW_Post_remove_header = RAW_Post(2:end,:);

% Format time and date of the data 
Time = extractAfter(RAW_Post_remove_header(:,3),"T");
Date = extractBefore(RAW_Post_remove_header(:,3),"T");

% Convert date format 
Date = datetime(Date,'InputFormat','yyyy-MM-dd');
formatOut = 'mm/dd/yyyy';
Date = string(datestr(Date,formatOut));
Date = cellstr(Date);
Time = string(datestr(Time, 'HH:MM:SS' ));
Time = cellstr(Time);

%%%% Addding new column for drug name 
% Drug_info = ["levetiracetam","lacosamide","lorazepam","phenytoin",...
%    "fosphenytoin","phenobarbital", "carbamazepine","valproate","divalproex","topiramate","clobazam",...
%    "lamotrigine", "oxcarbazepine","diazepam","gabapentin","dexmedetomidine","zonisamide","clonazepam","propofol","midazolam","ketamine","pentobarbital"];

% Drug data point of interest 
Drug_info = ["levetiracetam","lacosamide","lorazepam","phenytoin",...
   "fosphenytoin","phenobarbital", "carbamazepine","valproate","divalproex","topiramate","clobazam",...
   "lamotrigine", "oxcarbazepine","diazepam","zonisamide","clonazepam","propofol","midazolam","ketamine","pentobarbital"];


              



Drug_info = Drug_info';


%%%% Addding new column for drug name to save drug name 
drug_rename = cell(size(RAW_Post_remove_header,1),1);

%%%%% Convert drug units into upper case units for sorting purpose later

RAW_Post_remove_header(:,5) = cellfun(@upper,RAW_Post_remove_header(:,5),'UniformOutput',false);
% Find data belong to each drug 
 for ind_drug = 1: length(Drug_info)
     
        %%%% finding the string with indecies of the each unique drug and
        %%%% rename for better sorting in future 
     
     index =  ~cellfun('isempty',regexpi(RAW_Post_remove_header(:,10),Drug_info(ind_drug)));
     indx = find(index==1);
     drug_rename(indx,1) = cellstr(Drug_info(ind_drug));
 end  
 
 
 
% Save post epic data belong to each patient into their cell array  
for j = 1: size(Data_MNR_Sid,1)


indx = find(cell2mat(RAW_Post_remove_header(1:end,1))== cell2mat(Data_MNR_Sid(j,1)));
 column_index = [7 4 5];
 RAW_Buff{j,1} =[RAW_Post_remove_header(indx,1),Date(indx,1),Time(indx,1),RAW_Post_remove_header(indx, column_index),drug_rename(indx,1)];
    
end
 
% Look for pre epic data for each patient for each pre epic excel sheet 
for k =1:length(excelfiles )
    file = excelfiles(k);
    f = fullfile(file.folder,file.name);
    % read each pre epic excel file 
   [~,~,RAW]  =xlsread(f);
   
   %%% Converting time from number to string 
   str = string(datestr( cell2mat(RAW(2:end,6)), 'HH:MM:SS' ));
   RAW(2:end,6) = cellstr(str);
   
   RAW_local_Buff = RAW(2:end,:);
for j = 1: size(Data_MNR_Sid,1)
    % Comparing MNR so each patient data can be found save into particular
    % cell array 
   indx = find(cell2mat(RAW_local_Buff(1:end,1))== cell2mat(Data_MNR_Sid(j,1)));

   column_index = [1 5 6 12 13 14 15];
   % If data is found then concatenate cell array for the particular
   % patient
    RAW_Buff{j,1} =[ RAW_Buff{j,1};RAW_local_Buff(indx,column_index)];
end
end

% Check for patient whose data is not found 
Not_found = 1;
for check_ind = 1:size(RAW_Buff,1);
% Check if the patient found in the pre and post  epic
if isempty(RAW_Buff{check_ind,1})
    sprintf('%s_Not Found',Data_MNR_Sid{check_ind,3})
    Data_not_found(Not_found,:) = Data_MNR_Sid(check_ind,:);
    Not_found = Not_found+1;
else
    % Store data of the particular patient  into a excel sheet
      s = strsplit(Data_MNR_Sid{check_ind,3},{'sid'});
      file_Name = sprintf('sid%04s_MAR_Pre_Post_Epic.xlsx',s{2});
     %file_Name = sprintf('%s_MAR_Pre_Epic.xlsx',Data_MNR_Sid{check_ind,3});
     folder=['C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\Sid_Pre_Post_EPIC\'  file_Name];
     % Adding header to the data 
     Data_Header =[ RAW(1,1),RAW(1,5),RAW(1,6),RAW(1,12),RAW(1,13),RAW(1,14),RAW(1,15)];
     RAW_Buff_inc_header = [Data_Header ;RAW_Buff{check_ind,1}];
     xlswrite(folder, RAW_Buff_inc_header);
end
end


%%% Storing Data not found in excel sheet
file_Name = sprintf('Sids_Not_Found.xlsx');
folder=['C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Codes for Format CNN_Drug_Spike\Step_3_AED Extractions for 50 patients\'  file_Name];
   
Data_Not_Found = [Data_MRN_Header ;Data_not_found];
xlswrite(folder, Data_Not_Found);
 