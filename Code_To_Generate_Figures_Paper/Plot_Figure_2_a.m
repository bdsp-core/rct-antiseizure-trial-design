%%%% Calculating the seizure burden from the CNN labels %%%%
clc
clear all
close all
excelfiles = dir('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\RCT_sim_paper\Code\CNN_Label\*.csv')

%excelfiles = dir('C:\Users\CDAC_01\Dropbox (Partners HealthCare)\CausalModeling_IIIC\data_to_share\human_label_24h\*.csv')
% selected Drug for PK model
Drug_info = ["levetiracetam","lacosamide","midazolam","phenobarbital",...
    "pentobarbital","propofol","valproate"];
% half life
%  'lacosamide':[13], 'levetiracetam':[6], 'midazolam':[1.5], 'pentobarbital':[15]
% 'phenobarbital':[53],  'phenytoin':[22], 'propofol':[1.5],'valproate':[8]
t_half = [48,66,15,474,195,2,96]';
k = log(2)./t_half;

Drug_table = table(Drug_info',t_half ,k,'VariableNames',{'Drug_info','t_half','k'});


% check lognormal-drug_effect is non negative
y_log_drug_effect = {};

Sid = {};
%%% parameter buffer %%%
log_normal_parameters_buffer = [];
MA_parameters_buffer = [];
Length_seizure_burden = [];
%  for k =1:length(excelfiles )
fig = figure('Position', get(0, 'Screensize')) ;
cmap = colormap(lines(length(excelfiles )));
%  for k =1:length(excelfiles )
sid_check = 0;
for k =1:length(excelfiles )
    %k = [2,3,5,6,8,9,10,11];
    if (k~=12) && (k~=50)
        sid_check = sid_check+1;
        file = excelfiles(k);
        Sid{k}  = extractBefore(file.name,"_");
        f = fullfile(file.folder,file.name);
        RAW  =readtable(f,'ReadVariableNames',false);
        CNN_Label = table2array(RAW) ;% cell 2 matrix conversion
        Time_window = 10; % minutes
        No_min_Hr = 60; % Number of minutes in one hour
        Window_size = (Time_window*60)/2; % Window size is in seconds
        
        % Spike detection
        % Load spike mat file
        Temp_ = load(sprintf('C:\\Users\\CDAC_01\\Dropbox (Partners HealthCare)\\RCT_sim_paper\\Code\\Spike Data\\%s_Spike_Artifacts.mat',Sid{k}));
        
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
        %%% Finding seizures in the matrix
        IIC_logic = CNN_Label_window>=1 & CNN_Label_window <= 4 ;
        
        
        %%% Assign 1's where find seizure
        IIC_burden(IIC_logic) = 1;
        
        %%% Assign NaN's where their is NaN in the label data
        IIC_burden(logic_nan) = nan;
        
        %%% finding seizure burden by calculating the percentage of the mean at a
        %%% given time "seizure burden minutes / hour"
        
        IIC_final_burden = nanmean(IIC_burden,2) * 60;
        
        % End of IIC burden Processing
        %%% Calculate Drug Concentration %%%
        % load mat file of the sid
        Temp_ = load(sprintf('C:\\Users\\CDAC_01\\Dropbox (Partners HealthCare)\\RCT_sim_paper\\Code\\Combined_Drug_Normalize\\%s.mat',Sid{k}));

        
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
        
        
        testmat11 = py.importlib.import_module('Curve_fit_Python_ver16_Only_Log');
        %testmat12 = py.importlib.import_module('Curve_fit_Python_ver07_A_std_mu');
        %%%% x with respect to seizure burden
        
        x = (Time_window/No_min_Hr)*(1:length(IIC_final_burden));
        %%% Concatenating all NaN x for curve fitting
        % concatenate time and Drug Concatenation
        x_DC = [x',Drug_Conc'];
        x_curve_fit = x_DC(indices_Non_Nans,:);
        %%% Concatenating all NaN Seizure burden  for curve fitting
        IIC_final_burden_curve_fit = IIC_final_burden(indices_Non_Nans);
        
        Curvefit_params = testmat11.generate_Params(x_curve_fit,IIC_final_burden_curve_fit);
        
        
        
        %%% generate paramters from the python curve fit function
        parameters  = double(py.array.array('d',py.numpy.nditer(Curvefit_params)));
        
        %%%% Save those paramters for calculating the mean and covaraince for the
        %%%% three paramters mean, sigma , Peak Amplitude
        log_normal_parameters_buffer(k,:) = parameters ;
        
        
        %%% calcualate residuals %%%%
        y_final_burden_curve_fit = testmat11.f(x_curve_fit,parameters(1),parameters(2),parameters(3));
        
        
        
        % calculation of combining lognormal and Drug
        y1_log_Drug= testmat11.f(x_DC,parameters(1),parameters(2),parameters(3));
        y1_log_Drug= double(py.array.array('d',py.numpy.nditer(y1_log_Drug)));
        y_log_drug_effect{k} = y1_log_Drug;
        % calculating only lognormal
        y1_log= testmat11.f_log(x_DC(:,1),parameters(1),parameters(2),parameters(3));
        y1_log= double(py.array.array('d',py.numpy.nditer(y1_log)));
        
        y1_log(y1_log >60) = 60;
        y1_log(y1_log<0) = 0;
        
        
        % y1_noise_norm = y1_noise_norm *normalized_threshold;
        
        
        % %%% clipping negative amplitude to zero
        % y1_noise(y1_noise<0)=0;
        % %%% clipping  amplitude greater than 60 to 60
        % y1_noise(y1_noise>60)=60;
        
        Length_seizure_burden(k,1)= length(IIC_final_burden);
        %%%% Plot for highlighting NaN's values
        Y_NaN = NaN(size(IIC_final_burden));
        Y_NaN(indices ) = 65;
        
        
        
        
        
        
        
        
%         plot(x,y1_log,'LineWidth',1,'Color',cmap(k,:))
        plot(x,y1_log,'LineWidth',1,'Color',[0.9100    0.4100    0.1700])
        hold on
        
        
        % xlim([0 240])
        % ylim([-10 70])
        xlabel('Time in Hrs','FontSize',14)
        ylabel('Epileptiform Activity Bruden (min/h)','FontSize',14)
        %legend('SZ CNN Data ','Curve fit','NaN Val','IIC Gen from Log curve')
        % legend('IIC Burden ','Drug free IIC Burden','NaN Val')
        % title(sprintf('%s IIC burden and Drug free IIC Burden ', extractBefore(file.name,"_")))
        title('Epileptiform Activity Bruden','FontSize',14)
        
        
    end
end

saveas(fig,'C:\\Users\\CDAC_01\\Dropbox (Partners HealthCare)\\RCT_sim_paper\\Code\\Code_To_Generate_Figures_Paper\\Figure_2_a.jpeg');


% %%% Transpose Sid for data variable
%  Sid =Sid';
%  %%% Dictonary for looking at the details of the parameter details
%  %%% generated from curve fit function
% dictonary_Header = [{'Sids'},{'mu'},{'Sigma'},{'Peak_Value'},cellstr("DC_W" + (1:length(Drug_info))),{'Length Seizure Burden'},{'er_mean'},{'er_variance'},cellstr("lag_" + (1:no_lags))];
% dictonary_buffer = [Sid,num2cell(log_normal_parameters_buffer),num2cell(Length_seizure_burden),num2cell(MA_parameters_buffer)];
% dictonary_data = [dictonary_Header;dictonary_buffer];
% xlswrite('Parameters 50 patients without eps_DC_relu_Clip.xlsx',dictonary_data )
