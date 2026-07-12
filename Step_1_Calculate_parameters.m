%%%% Calculating the seizure burden from the CNN labels %%%%
%%%% Please run the this Step_1_Calculate_parameters code using anaconda
%%%% matlab
% step 1: first run anaconda prompt
% step 2 : Run matlab in anaconda prompt
% step 3: run Step_1_Calculate_parameters code using matlab open through
% anaconda prompt
clc
clear all
close all
excelfiles = dir('CNN_Label\*.csv')

% selected Drug for PK model
Drug_info = ["levetiracetam","lacosamide","midazolam","phenobarbital",...
    "pentobarbital","propofol","valproate"];
% half life with respect to minutes
% 'lacosamide':[66],     #  11h, (5-15h)
%         'levetiracetam':[48],  #  8h
%         'midazolam':[15],      #  2.5h
%         'pentobarbital':[195], # 32.5h (15-50h)
%         'phenobarbital':[474], # 79h
%         'phenytoin':[147],     # 24.5h (7-42h)
%         'propofol':[2],        # 20minutes (3-12h after long time) (needs 3 differential equations)
%         'valproate':[96]       # 16h

t_half = [48,66,15,474,195,2,96]';
k = log(2)./t_half;

Drug_table = table(Drug_info',t_half ,k,'VariableNames',{'Drug_info','t_half','k'});


% check lognormal-drug_effect is non negative
y_log_drug_effect = {};

Sid = {};
Sid_buff ={};
%%% parameter buffer %%%
log_normal_parameters_buffer = [];

Length_seizure_burden = [];

%%% variable to store cdf
cdf_SB = [];

Sid_store = 1;
for k =1:length(excelfiles )
    
    if (k~=12) && (k~=50)  % sid0716 & sid1881 are not selected due to high standard deviation
        
        file = excelfiles(k);
        Sid{k}  = extractBefore(file.name,"_");
        Sid_buff{Sid_store}  = extractBefore(file.name,"_");
        f = fullfile(file.folder,file.name);
        RAW  =readtable(f,'ReadVariableNames',false);
        CNN_Label = table2array(RAW);% cell 2 matrix conversion
        Time_window = 10; % minutes
        No_min_Hr = 60; % Number of minutes in one hour
        Window_size = (Time_window*60)/2; % Window size is in seconds
        
        % Spike detection
        % Load spike mat file
        Temp_ = load(sprintf('Spike Data\\%s_Spike_Artifacts.mat',Sid{k}));
        Artifacts = Temp_.Buff_art;
        Spikes_ = Temp_.Buff_yp;
        % Convert Spikes with Artifacts to NaN's
        Spikes_(Artifacts ==1)=nan;
        Spikes_Window = [];
        index = 1;
        for i = 1 : Window_size : length(Spikes_)-(Window_size-1)
            Spikes_Window(index,:) =Spikes_(i:i+(Window_size-1));
            index = index + 1;
        end
        Spike_rate = nanmean(Spikes_Window,2);
        
        % End of Spike processing
        
        % Begin of IIC burden Processing
        CNN_Label(Artifacts ==1)=nan;
        CNN_Label_window = [];
        index = 1;
        for i = 1 : Window_size : length(CNN_Label)-(Window_size-1)
            CNN_Label_window(index,:) =CNN_Label(i:i+(Window_size-1));
            index = index + 1;
        end
        
        %%% Find all nan in CNN_labels_window to replace it to sz_burden_matrix
        logic_nan =  isnan(CNN_Label_window);
        %%%% Sezure matrix of zeros of size eequal to CNN Label matrix
        IIC_burden = zeros(size(CNN_Label_window));
        %%% Finding seizure,LPD,LRPDA,GPD in the matrix
        IIC_logic = CNN_Label_window>=1 & CNN_Label_window <= 4 ;
        
        
        %%% Assign 1's where find seizure,LPD,LRPDA,GPD
        IIC_burden(IIC_logic) = 1;
        
        %%% Assign NaN's where their is NaN in the label data
        IIC_burden(logic_nan) = nan;
        
        %%% finding seizure burden by calculating the percentage of the mean at a
        %%% given time "seizure burden minutes / hour"
        
        IIC_final_burden = nanmean(IIC_burden,2) * 60;
        
        % End of IIC burden Processing
        
        %%% Calculate Drug Concentration %%%
        % load Drug Concentration mat file of the sid
        Temp_ = load(sprintf('Combined_Drug_Normalize\\%s.mat',Sid{k}));
        
        Drug_name = Temp_.Drug_info;
        Drug_Buffer = Temp_.Drugs_Normalize;
        Drug_window = [];
        index = 1;
        for i = 1 : Window_size : length(Drug_Buffer)-(Window_size-1)
            Drug_window(index,:) =mean(Drug_Buffer((i:i+(Window_size-1)),:));
            index = index + 1;
        end
        
        % PK model for the
        
        % selct the drugs we are applying PK model
        [~,pos]=ismember(Drug_info,Drug_name);
        Drug_=Drug_window(:,pos);
        % exponential of drugs over time
        exp_k = exp(-Drug_table.k .*(0:size(Drug_,1)-1));
        % Drug concentration over a time
        Drug_Conc = [];
        for drug_id = 1:size(Drug_,2)
            % performing convolution of each drug with respect to elimation rate
            Drug_Conc(drug_id,:)=conv(Drug_(:,drug_id),exp_k(drug_id,:));
        end
        Drug_Conc = Drug_Conc(:,1:size(Drug_,1));
        
        
        %%%% Looking for NaN's in the seizure
        indices = find(isnan(IIC_final_burden) == 1);
        indices_Non_Nans = find(~isnan(IIC_final_burden) == 1);
        %%%% Load python curve fitting function in matlab
        
        %%% Load curve fit code that is written in python language.
        testmat = py.importlib.import_module('Curve_fit_Python');
        
        %%%% x with respect to seizure burden
        
        x = (Time_window/No_min_Hr)*(1:length(IIC_final_burden));
        %%% Concatenating all non NaN x for curve fitting
        % concatenate time and Drug Concatenation
        x_DC = [x',Drug_Conc'];
        x_curve_fit = x_DC(indices_Non_Nans,:);
        %%% Concatenating all non NaN IIC burden  for curve fitting
        IIC_final_burden_curve_fit = IIC_final_burden(indices_Non_Nans);
        
        Curvefit_params = testmat.generate_Params(x_curve_fit,IIC_final_burden_curve_fit);
        
        
        
        %%% generate paramters from the python curve fit function
        parameters  = double(py.array.array('d',py.numpy.nditer(Curvefit_params)));
        
        %%%% Save those paramters for calculating the mean and covaraince for the
        %%%% three paramters mean, sigma , Peak Amplitude
        log_normal_parameters_buffer(Sid_store,:) = parameters ;
        
        
        % calculation of combining lognormal and Drug
        y1_log_Drug= testmat.f(x_DC,parameters(1),parameters(2),parameters(3),parameters(4),parameters(5),parameters(6),parameters(7),parameters(8),parameters(9),parameters(10));
        y1_log_Drug= double(py.array.array('d',py.numpy.nditer(y1_log_Drug)));
        y_log_drug_effect{Sid_store} = y1_log_Drug;
        
        % calculating only lognormal
        y1_log= testmat.f_log(x_DC(:,1),parameters(1),parameters(2),parameters(3));
        y1_log= double(py.array.array('d',py.numpy.nditer(y1_log)));
        
        y1_log(y1_log >60) = 60;
        y1_log(y1_log<0) = 0;
        
        
        %%%% finding the conditiional probability of the data less than 14 days.
        cdf_SB(Sid_store,1)=logncdf(24*14,parameters(1),parameters(2));
        
        Length_seizure_burden(Sid_store,1)= length(IIC_final_burden);
        %%%% Plot for highlighting NaN's values
        Y_NaN = NaN(size(IIC_final_burden));
        Y_NaN(indices ) = 65;
        
        
        
        fig = figure('Position', get(0, 'Screensize'))
        %figure
        subplot(4,1,1);
        plot(x,IIC_final_burden,'LineWidth',1)
        hold on
        plot(x,y1_log_Drug,'LineWidth',1)
        hold on
        plot(x,Y_NaN,'LineWidth',5,'color',[0 0 0])
        ylim([-10 70])
        xlabel('Time in Hrs')
        ylabel('seizure burden (min/h)')
        
        legend('SZ CNN Data ','lognormal-DrugCon','NaN Val')
        title(sprintf('%s Seizure Burden Lognormal Calculation ', extractBefore(file.name,"_")))
        
        subplot(4,1,2);
        colorstring = 'rmckgby';
        loc_drugs = find(sum(Drug_Conc,2)>0);
        for i = 1:length(loc_drugs)
            plot(x,Drug_Conc(loc_drugs(i),:),colorstring(loc_drugs(i)),'LineWidth',1)
            hold on
        end
        hold off
        %ylim([-10 70])
        xlabel('Time in Hrs')
        ylabel('Normalized Drug Concentration')
        legend(Drug_info (loc_drugs))
        title(sprintf('%s Drug Concentration ', extractBefore(file.name,"_")))
        
        subplot(4,1,3);
        plot(x,IIC_final_burden,'LineWidth',1)
        hold on
        plot(x,y1_log,'LineWidth',1)
        hold on
        plot(x,Y_NaN,'LineWidth',5,'color',[0 0 0])
        ylim([-10 70])
        xlabel('Time in Hrs')
        ylabel('seizure burden (min/h)')
        legend('SZ CNN Data ','lognormal','NaN Val')
        title(sprintf('%s Seizure Burden Lognormal Calculation ', extractBefore(file.name,"_")))
        
        subplot(4,1,4);
        plot(x,Spike_rate,'LineWidth',1)
        ylim([0 1])
        xlabel('Time in Hrs')
        ylabel('Spike Rate')
        legend('Spike Rate')
        title(sprintf('%s Spike Rate ', extractBefore(file.name,"_")))
        
        saveas(fig,sprintf('.\\Sim_fit_figures\\%s IIC_Burden_fit.jpeg', extractBefore(file.name,"_")))
        Sid_store = Sid_store+1;
    end
end

% %%% Transpose Sid for data variable
Sid =Sid';
Sid_buff = Sid_buff';
%%% Dictonary for looking at the details of the parameter details
%%% generated from curve fit function
dictonary_Header = [{'Sids'},{'mu'},{'Sigma'},{'Peak_Value'},cellstr("DC_W" + (1:length(Drug_info))),{'Length Seizure Burden'},{'CDF<14days'}];
dictonary_buffer = [Sid_buff,num2cell(log_normal_parameters_buffer),num2cell(Length_seizure_burden),num2cell(cdf_SB)];
dictonary_data = [dictonary_Header;dictonary_buffer];
xlswrite('Parameters 48 patients final.xlsx',dictonary_data )
