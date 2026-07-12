%%% Data Analysis MAR Pre Epic Data %%%%

clc
clear all
close all

% Find all 50 sid names thatr are included in the study


Sid = ["sid0742","sid0861","sid0877","sid0955","sid0970","sid0351","sid0959","sid0920","sid0345","sid1127",...
    "sid1125","sid0868","sid1158","sid0716","sid0708","sid0699","sid0341","sid1535","sid1881","sid0353","sid1144",...
    "sid1577","sid1532","sid1486","sid1832","sid0834","sid0854","sid0855","sid1211","sid1118","sid0350","sid1194",...
    "sid1119","sid0846","sid0871","sid0320","sid0949","sid0804","sid1109","sid0070","sid1112","sid0962","sid0287",...
    "sid0286","sid0967","sid1117","sid0826","sid0922","sid0958","sid0850"];


%%%% Load excel sheet for selcting the CNN data between the admit and
%%%% Discharge dates of the patient

[RAW_sid_SAGE_SAH]  =readtable('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\AED_TRIAL_PROJECT\BadChannelDetetction_demo\BadChannelDetetction_demo\SAH\Info 65 clients in SAH SAGE DataBase.xlsx');
RAW_sid_SAGE_SAH = table2cell(RAW_sid_SAGE_SAH);

% loop for processing data for each patient

for j = 1:length(Sid)
    
    % buffer to store CNN data
    Buff = [];
    
    % Load xlsx files of CNN label for the particular SID
    
    excelfiles = dir(sprintf('Z:\\Projects\\Rajesh\\mgh_iiica_online\\SAH\\%s\\%s_Labels\\*.csv',Sid(j),Sid(j)));
    start_time_Human_Buff = {};
    RAW_Buff_Buff = {};
    Sid(j)
    %%%% Saving first CNN label Date_time  for saving the file name
    save_file_name = excelfiles(1);
    save_file_name  = save_file_name.name;
    save_file_name  = extractBetween(save_file_name  ,"online_",".csv");
    save_file_name = extractAfter(save_file_name,"_");
    
    % Each patient has multiple cEEG data for each encounter
    % loop for saving CNN data and extract of start time of their CNN data
    
    for k =1:length(excelfiles )
        file = excelfiles(k);
        
        f = fullfile(file.folder,file.name);
        [~,~,RAW]  =xlsread(f);
        
        % Extracting start date and time of each segment of CNN data for cEEG
        EEG_Start_Time = extractAfter(file.name,"_");
        EEG_Start_Time = extractAfter(EEG_Start_Time ,"_");
        start_time = extractBefore(EEG_Start_Time ,"_Bad");
        % Convert into datetime format and Save CNN label start time for each segment
        start_time_Human_Buff{k} = datetime(start_time,'InputFormat','yyyyMMdd_HHmmss');
        %%% Insert NaNs starting and Ending of the CNN label due to CNN classfier did not
        % include first and last two seconds
        
        RAW_Buff_Buff{k} = [nan(2,1); cell2mat(RAW);nan(2,1)];
        
        
        
    end
    
    %%%% Select CNN data which is within Hospital admit and discharge dates
    
    % Find patient Hospital admit and discharge infomration
    indx = find(strcmp(RAW_sid_SAGE_SAH(:,3),sprintf('sid%04d',str2double(extractAfter(Sid(j),"sid")))))
    Hospital_admit =datenum( datetime(RAW_sid_SAGE_SAH{indx,5},'InputFormat','MM/dd/yyyy'));
    Hospital_discharge =datenum( (datetime(RAW_sid_SAGE_SAH{indx,6},'InputFormat','MM/dd/yyyy')+days(1)));
    
    %%% Convert start time into datetime format
    Start_time_EEG_datenum = cellfun(@(s) datenum(s), start_time_Human_Buff, 'UniformOutput', false);
    % convert cell to mat
    Start_time_EEG_datenum = cell2mat( Start_time_EEG_datenum);
    
    % Find CNN data which is within hospital admit and discharge dates
    indx_EEG_match = find( Start_time_EEG_datenum(:) >= Hospital_admit &  Start_time_EEG_datenum(:)< Hospital_discharge);
    start_time_Human = start_time_Human_Buff(indx_EEG_match);
    % Save CNN data which is within hospital admit and discharge dates.
    RAW_Buff  = RAW_Buff_Buff(indx_EEG_match);
    
    % concatenate the segments of CNN labelled data into one array and if there is a gap
    % between the segments that is filled with NaN values.
    
    for k1 = 1:length( start_time_Human)
        
        % Check for last segment so no need to apeend NaN values
        if k1 == length( start_time_Human)
            % Buffer to store concatenate CNN data
            Buff = [Buff,RAW_Buff{k1}'];
        else
            % Find end date for each segment
            end_time = start_time_Human{k1} + seconds(size(RAW_Buff{k1},1)*2);
            % start date time of next segment is greater than end date time of
            % present segment then append NaNs back of the present segment.
            if seconds(start_time_Human{k1+1} - end_time) > 0
                null_time_steps = round(seconds(start_time_Human{k1+1} - end_time)/2)
                Data_Time_step = [RAW_Buff{k1}',nan* ones(1,null_time_steps)];
                Buff = [Buff,Data_Time_step];
            else
                % if there is no gap between end date time of the present
                % segment and start date time of the next segment then no need
                % to append NaN values.
                Buff = [Buff,RAW_Buff{k1}'];
            end
        end
        
        
    end
    % Save continuous CNN labelled data which is concatenated by filling the missing data with NaNs
    folder_name = [ 'C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\CNN_Label\' ];
    filename = sprintf('%s_%s_CNN_Label.csv',Sid(j), save_file_name{1});
    folder=[folder_name  filename];
    csvwrite(folder, Buff');
end
