%%% Data Analysis MAR Pre Epic Data %%%%

clc
clear all
close all

% Load spike data calculated from cEEG using Spike detector
files_Sid = dir('Z:\Projects\Rajesh\Spike Detection\SSD_demo\SSD\spikewise\sid*.*')

% sampling frequency of spike
Fs_ = 128;

%%%% Load excel sheet for selcting the Spike data between the admit and
%%%% Discharge dates of the patient

[RAW_sid_SAGE_SAH]  =readtable('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\AED_TRIAL_PROJECT\BadChannelDetetction_demo\BadChannelDetetction_demo\SAH\Info 65 clients in SAH SAGE DataBase.xlsx');
RAW_sid_SAGE_SAH = table2cell(RAW_sid_SAGE_SAH);

% loop for processing data for each patient

for j = 1:length(files_Sid)
    
    % buffer to store spike detection
    Buff_yp = [];
    % buffer to store artifacts
    Buff_art = [];
    
    % Load xlsx files of Spike data for the particular SID
    
    matfiles = dir(sprintf('Z:\\Projects\\Rajesh\\Spike Detection\\SSD_demo\\SSD\\spikewise\\%s\\*.mat',files_Sid(j).name));
    Spike_start_time = {};
    RAW_Buff_Buff = {};
    files_Sid(j).name
    %%%% Saving first spike label Date_time  for saving the file name
    save_file_name = matfiles(1);
    save_file_name  = save_file_name.name;
    save_file_name  = extractBetween(save_file_name  ,"SSD_",".mat");
    %save_file_name = extractAfter(save_file_name,"_");
    
    % Each patient has multiple cEEG data for each encounter
    % loop for saving spike data and extract of start time of their spike data
    for k =1:length(matfiles )
        file = matfiles(k);
        data = load(fullfile(file.folder,file.name));
        % Extracting start date and time of each segment of spikes for cEEG
        EEG_Start_Time = extractAfter(file.name,"_");
        EEG_Start_Time = extractAfter(EEG_Start_Time,"_");
        EEG_Start_Time = extractBefore(EEG_Start_Time ,".mat");
        % Convert into datetime format and Save spike detector start time for each segment
        Spike_start_time{k} = datetime(EEG_Start_Time,'InputFormat','yyyyMMdd_HHmmss');
        % Save spike data for segment in the Buffer
        RAW_Buff_Buff{k} = data ;
        
        
        
    end
    
    %%%% Select spike data which is within Hospital admit and discharge dates
    % Find patient Hospital admit and discharge infomration
    indx = find(strcmp(RAW_sid_SAGE_SAH(:,3),files_Sid(j).name));
    Hospital_admit =datenum( datetime(RAW_sid_SAGE_SAH{indx,5},'InputFormat','MM/dd/yyyy'));
    Hospital_discharge =datenum( (datetime(RAW_sid_SAGE_SAH{indx,6},'InputFormat','MM/dd/yyyy')+days(1)));
    
    %%% Convert start time into datetime format
    Start_time_EEG_datenum = cellfun(@(s) datenum(s), Spike_start_time, 'UniformOutput', false);
    % convert cell to mat
    Start_time_EEG_datenum = cell2mat( Start_time_EEG_datenum);
    
    % Find spike data which is within hospital admit and discharge dates
    indx_EEG_match = find( Start_time_EEG_datenum(:) >= Hospital_admit &  Start_time_EEG_datenum(:)< Hospital_discharge);
    start_time_Spike = Spike_start_time(indx_EEG_match);
    
    % Save spike data which is within hospital admit and discharge dates.
    RAW_Buff  = RAW_Buff_Buff(indx_EEG_match);
    yp_2Sec =  {};
    art_2Sec =  {};
    % converting spike detection and artifacts into 2 seconds windows
    Window_size = 2*Fs_; % Window size is in seconds ie 2 seconds window
    % Loop for processing each cEEG segment of spike data
    for i1= 1:length(RAW_Buff)
        yp_window =[];
        art_window = [];
        index = 1;
        for i = 1 : Window_size : length(RAW_Buff{1,i1}.yp)-(Window_size-1)
            yp_window(index,1) =nanmax(RAW_Buff{1,i1}.yp(1,(i:i+(Window_size-1))));
            art_window(index,1) =nanmax(RAW_Buff{1,i1}.artifact(1,(i:i+(Window_size-1))));
            index = index + 1;
        end
        yp_window(yp_window>=0.43)=1;
        yp_window(yp_window<0.43)=0;
        yp_2Sec{i1}=yp_window;
        art_2Sec{i1}=art_window;
        
        
        
    end
    
    % Once we convert the spike data into 2 second spike data now need to
    % concatenate the segments of spike data into one array and if there is a gap
    % between the segments that is filled with NaN values.
    
    for k1 = 1:length( start_time_Spike)
         % Check for last segment so no need to apeend NaN values 
        if k1 == length( start_time_Spike)
            
            % Buffer to store spike info
            Buff_yp = [Buff_yp,yp_2Sec{k1}'];
            % Buffer to store artifacts info
            Buff_art = [Buff_art,art_2Sec{k1}'];
        else
            % Find end date for each segment
            end_time = start_time_Spike{k1} + seconds(size(yp_2Sec{k1},1)*2);
            % start date time of next segment is greater than end date time of
            % present segment then append NaNs back of the present segment.
            if seconds(start_time_Spike{k1+1} - end_time) > 0
                null_time_steps = round(seconds(start_time_Spike{k1+1} - end_time)/2)
                Data_Time_step_yp = [yp_2Sec{k1}',nan* ones(1,null_time_steps)];
                Data_Time_step_art = [art_2Sec{k1}',nan* ones(1,null_time_steps)];
                Buff_yp = [Buff_yp,Data_Time_step_yp];
                Buff_art = [Buff_art,Data_Time_step_art];
            else
                % if there is no gap between end date time of the present
                % segment and start date time of the next segment then no need
                % to append NaN values.
                Buff_yp = [Buff_yp,yp_2Sec{k1}'];
                Buff_art = [Buff_art,art_2Sec{k1}'];
            end
        end
        
        
    end
    % Savecontinuous spike data which is converted into 2 seconds segments and
    % concatenated by filling the missing data with NaNs
    save_path = [ 'C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\Spike Data\' ];
    filename = sprintf('%s_Spike_Artifacts.mat',files_Sid(j).name);
    
    save([save_path, filename],...
        'Buff_yp', 'Buff_art')
end
