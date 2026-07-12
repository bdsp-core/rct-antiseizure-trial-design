%%% Data Analysis MAR Pre Epic Data %%%%

clc
clear all
close all

Buff = cell(1,1);
Sid = ["sid0742","sid0861","sid0877","sid0955","sid0970","sid0351","sid0959","sid0920","sid0345","sid1127",...
    "sid1125","sid0868","sid1158","sid0708","sid0699","sid0341","sid1535","sid0353","sid1144",...
    "sid1577","sid1532","sid1486","sid1832","sid0834","sid0854","sid0855","sid1211","sid1118","sid0350","sid1194",...
    "sid1119","sid0846","sid0871","sid0320","sid0949","sid0804","sid1109","sid0070","sid1112","sid0962","sid0287",...
    "sid0286","sid0967","sid1117","sid0826","sid0922","sid0958","sid0850"];

fig = figure('Position', get(0, 'Screensize')) ;
Buff = cell(size(Sid,2),1);

[RAW_sid_SAGE_SAH]  =readtable('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\AED_TRIAL_PROJECT\BadChannelDetetction_demo\BadChannelDetetction_demo\SAH\Info 65 clients in SAH SAGE DataBase.xlsx');
 RAW_sid_SAGE_SAH = table2cell(RAW_sid_SAGE_SAH);

for j =1:length(Sid)

% Load xlsx files MAR sheet

excelfiles = dir(sprintf('Z:\\Projects\\Rajesh\\mgh_iiica_online\\SAH\\%s\\%s_Labels\\*.csv',Sid(j),Sid(j)));
start_time_Human_Buff = {};
RAW_Buff_Buff = {};
Sid(j)
for k =1:length(excelfiles )
    file = excelfiles(k);
    f = fullfile(file.folder,file.name);
   [~,~,RAW]  =xlsread(f);
   
   file_Name = extractAfter(file.name,"_");
   file_Name = extractAfter(file_Name ,"_");
   start_time = extractBefore(file_Name ,"_Bad");
   start_time_Human_Buff{k} = datetime(start_time,'InputFormat','yyyyMMdd_HHmmss');
   RAW_Buff_Buff{k} = [6*ones(2,1); cell2mat(RAW);6*ones(2,1)];
  
end
%%%% Select start_time between the admit and discharge dates 
 indx = find(strcmp(RAW_sid_SAGE_SAH(:,3),sprintf('sid%04d',str2double(extractAfter(Sid(j),"sid")))))
 Hospital_admit =datenum( datetime(RAW_sid_SAGE_SAH{indx,5},'InputFormat','MM/dd/yyyy'));
 Hospital_discharge =datenum( (datetime(RAW_sid_SAGE_SAH{indx,6},'InputFormat','MM/dd/yyyy')+days(1)));
 
 %%% Convert start time into datetime format 
Start_time_EEG_datenum = cellfun(@(s) datenum(s), start_time_Human_Buff, 'UniformOutput', false);
 % convert cell to mat 
Start_time_EEG_datenum = cell2mat( Start_time_EEG_datenum);
indx_EEG_match = find( Start_time_EEG_datenum(:) >= Hospital_admit &  Start_time_EEG_datenum(:)< Hospital_discharge);
start_time_Human = start_time_Human_Buff(indx_EEG_match);
RAW_Buff  = RAW_Buff_Buff(indx_EEG_match);


%%%%% May be this block need to modify based on the chnage of the code in

for k1 = 1:length( start_time_Human)
    if k1 == length( start_time_Human)
         Buff{j,1} = [Buff{j,1},RAW_Buff{k1}'];
    else
        end_time = start_time_Human{k1} + seconds(size(RAW_Buff{k1},1)*2);
        if seconds(start_time_Human{k1+1} - end_time) > 0
         null_time_steps = round(seconds(start_time_Human{k1+1} - end_time)/2)
         Data_Time_step = [RAW_Buff{k1}',6* ones(1,null_time_steps)];
         Buff{j,1} = [Buff{j,1},Data_Time_step];
        else 
          Buff{j,1} = [Buff{j,1},RAW_Buff{k1}'];   
        end
    end
        
    
end
end

%%% Finding the maxmium length sid 
max_length  = 0;
for i = 1:size(Buff,1)
max_length = max([size(Buff{i,1},2),max_length]);
end
Buffer_converted = [];

%%% Making sure all sids have same length for imagesc function if of
%%% different length then append zeros to rest of the sid when compared to
%%% max length sid
for i = 1:size(Buff,1)
  %% GAP is going with 6 label
    Buffer_converted(i,:) = [ Buff{i,1},6*ones(1,max_length-size(Buff{i,1},2))];
    
end

% Plotting swimmer plot 
% figure
cmap = zeros(7, 3);

cmap = [0, 0, 1; ...       % other - dark blue for 0
  1, 0, 0; ...   % seizure - red for 1
  1, 0.5843, 0; ...       % lpd - orange for 2
  1, 1, 0; ...   % gpd - yellow for 3
  0, 1, 0; ...       % lrda - green for 4
  0.5843    0.8157    0.9882; ...       % grda - light blue for 5
  0.7, 0.7, 0.7]      % gap - grey for NaN NaN was label as 6
 


% Display the image with the colormap applied.
imagesc(Buffer_converted)
% Zoom in the graph with respecgt to time
%%%% 1 no zoom if scale increases zoom increases 1----> 10 [1 no zoom 10
%%%% resolution /10]
Zoom_res = 1;


xticklabels = (2*(linspace(0, max_length/Zoom_res, 10)))/3600;
xticks = linspace(1, max_length/Zoom_res, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
set(gca, 'Ytick',1:size(Buffer_converted,1),'YTickLabel',Sid);

xlim([0  max_length/Zoom_res])
colormap(cmap)
colorbar('Ticks',[0,1,2,3,4,5,6],...
         'TickLabels',{'Other','Seizure','LPD','GPD','LRDA','GRDA','GAP'})
xlabel('Time in hours')
ylabel('patient')   
saveas(fig,'C:\\Users\\CDAC_01\\Dropbox (Partners HealthCare)\\AED_TRIAL_PROJECT\\BadChannelDetetction_demo\\BadChannelDetetction_demo\\SAH\\50 Clients Point of interest\\Testing simulation with out EPS\\Code_Figures_Paper\\Swimmer_Plot.jpeg');