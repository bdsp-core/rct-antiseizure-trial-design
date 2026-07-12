%%%% Calculating the seizure burden from the CNN labels %%%%
clc
clear all
close all
%% Input to claculate sample size and effect size

T_delay_Array = [1,14,25,40,50,60,70]; % time delay drug administration in hours

Monitoring_Offset = 72 ;% 24 hrs   ; % Monitoring offset time in hours

% Selecting the ratio of patient in treatment group with Complete cessation value between  between 0 and 1
% 60% means 0.6, 50% means 0.50 and 75% means 0.75
treatment_group_pop_sel = 0.60;

% delay in drug response
T_drug_res = 14; % time in hours example [0,1,14,24] hours

% Selecting outcome
% Out comes "tSB","pSB12", "pPSB","pSB24"

%%%% Seizure burden means EA burden %%%%

% tSB total seizure burden from drug intervention
% pSB72 Post seizure burden for 72 hours from drug intervention
% pPSB Post seizure burden for 24 hours from peak after drug intervention
% pSB24 Post seizure burden for 24 hours from drug intervention
% Out_Come_array = ["tSB","pSB72", "pPSB","pSB24"];

Out_Come_array = ["tSB","pSB72","pSB24"];

Seizure_time_course_days = 14; % length of seizure time course
Time_window = 10; % minutes
No_min_Hr = 60; % Number of minutes in one hour

Log_fit_para = readtable('Parameters 48 patients final.xlsx');

%%
%%%% save only the mean , standard deviation and amplitude parameters of the lognormal into a buffer

parameters_buffer=[Log_fit_para.mu,Log_fit_para.Sigma,Log_fit_para.Peak_Value];


%%% Convert Std into positive due to std^2 in lognormal curvefit may pick
%%% negative values
parameters_buffer(:,2)=abs(parameters_buffer(:,2));
% Parameters are saved into parameters_buffer_v1
parameters_buffer_v1 = parameters_buffer;




% Initial sample and Effective sizes array

% Sample buffer to save calculated sample size for median, lower 2.5
% quartile and upper 2.5 quartile

Sample_size_Median = [];
Sample_size_Low = [];
Sample_size_Upper =[];

% Sample buffer to save calculated rounded sample size to nearest even number for median, lower 2.5
% quartile and upper 2.5 quartile
Sample_size_Median_round=[];
Sample_size_Low_round = [];
Sample_size_Upper_round =[];

% Buffer to save effect sizes for different outcomes
Effect_Size_Median = [];
Effect_Size_Low  = [];
Effect_Size_Upper =[];

% For loop for each outcomes
for i_outcome = 1:length(Out_Come_array)
    Out_Come = Out_Come_array(i_outcome); % Buffer to save each outcome
    
    % inner for loop  for each treatment delay
    for itr_time = 1:length(T_delay_Array)
        % Buffer for saving sample size and effect size for every outcome
        % at every time delay
        Sample_Size_Buff = [];
        Effective_Size_Buff = [];
        % calcuate drug efficacy
        % time of delay of treatment
        T_delay = T_delay_Array(itr_time);% Buffer to save each treatment delay
        
        Drug_Effect_Placebo = 0; % No Attenuation for Placebo
        
        ite_n = 1; % initialising boost variable
        check_it = 1 ; % check for test
        
        % number of iterations to run boostrapping step default is 1000
        %         while ite_n <=1000
        while ite_n <=1000
            % input paramters
            check_gen = 1;
            
            Select_no_Trails = 800000;
            % while loop to check total simulations is greater than 100000
            while check_gen == 1
                check_gen = check_gen+1;
                
                %%% calculation of mean and covariance of the data
                
                % boost strapping the data
                parameters_buffer_v2 = parameters_buffer_v1(randsample(size(parameters_buffer_v1,1),size(parameters_buffer_v1,1),true),1:3);
                para_mean = mean(parameters_buffer_v2);
                C = cov(parameters_buffer_v2);
                %%% Randomly generating the n sample using Multivariate normal random numbers
                
                Generate_burden_parameters = mvnrnd(para_mean,C,Select_no_Trails);
                
                
                %%%%% Select mu and sigma based within boostrapped data
                
                mu_select = find(Generate_burden_parameters(:,1)>=min(parameters_buffer_v2(:,1)) & Generate_burden_parameters(:,1)<=max(parameters_buffer_v2(:,1)));
                Generate_burden_parameters = Generate_burden_parameters (mu_select,:);
                % Sigma of select based on within boostrapped data
                sigma_select = find(Generate_burden_parameters(:,2)>=min(parameters_buffer_v2(:,2)) & Generate_burden_parameters(:,2)<=max(parameters_buffer_v2(:,2)));
                Generate_burden_parameters = Generate_burden_parameters (sigma_select ,:);
                
                
                % select the patients whose cummulative distribution
                % function is 90% percent within 14 days can also ratio
                cdf_14days = find(logncdf((24*Seizure_time_course_days),Generate_burden_parameters(:,1),Generate_burden_parameters(:,2))>0.90);
                Generate_burden_parameters = Generate_burden_parameters (cdf_14days,:);
                %%% Select peaks which are in interval in 0 to 60
                peak_select = find(Generate_burden_parameters(:,3)>=0 & Generate_burden_parameters(:,3)<=60);
                Generate_burden_parameters = Generate_burden_parameters (peak_select,:);
                
                % X axis seizure burden taxis for each 10 minutes
                Seizure_burden_time_course = (Time_window/No_min_Hr):(Time_window/No_min_Hr):24*Seizure_time_course_days;
                
                % lognormal curve
                % Repeat matrix with seizure burden time courses for
                % parallel calculating the seizure burden
                Seizure_burden_time_course_repeat = repmat(Seizure_burden_time_course,size(Generate_burden_parameters,1),1);
                
                
                fun_curve  = Generate_burden_parameters(:,3)'.*exp(-(((log(Seizure_burden_time_course_repeat')-Generate_burden_parameters(:,1)')./(sqrt(2)*Generate_burden_parameters(:,2)')).^2));
                
                % remove negative values from the data
                fun_curve = fun_curve(:,find(sum(fun_curve,1)>0));
                % just a cross check to select simulations greater than 0
                fun_curve = fun_curve(:,find(min(fun_curve, [], 1)>=0));
                % just a cross check to select max simulation value equal to or
                % less than 60
                fun_curve = fun_curve(:,find(max(fun_curve, [], 1)<=60));
                
                % find seizure onset time and seizure on time + monitor window is
                % less than total seizure time course.
                seizure_onset_group = [];
                for group_ind = 1:size(fun_curve,2)
                    seizure_onset_group(group_ind,1) = find(fun_curve(:, group_ind)>0,1);
                end
                % Select seizure burdens atleast 96 hours from seizure
                % onset.
                fun_curve = fun_curve(:,find((seizure_onset_group+((24+Monitoring_Offset)*(No_min_Hr/Time_window)))<size(fun_curve,1)));
                
                % Look for patients with atleast time duration of 30 hours from peak
                [peak_Amp,peak_Pos]=max(fun_curve, [], 1);
                fun_curve = fun_curve(:,find((peak_Pos+((No_min_Hr/Time_window)*30))<size(fun_curve,1)));
                
                generate_Y_buffer = fun_curve';
                
                generate_simulations = generate_Y_buffer(find(sum(generate_Y_buffer,2)~=0),:);
                % If condition to check total number of simulations greater
                % than 100000 if less than while loop runs again
                if size(generate_simulations,1)<100000
                    check_gen  = 1;
                    Select_no_Trails = Select_no_Trails+round(0.25*Select_no_Trails); % Increase total number of trails
                end
            end % end of while loop for generating simulations
            
            
            check_it = check_it+1 ;
            
            % Suffle the data so randomly select 100000 simulations
            Shuffle_indices_group = randperm(size(generate_simulations,1));
            
            
            % Random select 100000 simulations  so each group has 50,000 simulations
            generate_Y_noise_buffer_RCT = generate_simulations (Shuffle_indices_group(1:100000),:);
            
            % Shuffle the data so can split into placebo and intervension/First line AED group
            Shuffle_indices = randperm(size(generate_Y_noise_buffer_RCT,1));
            % randomly selecting the 50% percentage between Intervension and placebo
            % group
            idx_sub_1st_Line=ismember(Shuffle_indices,Shuffle_indices(1:round(0.50*size(generate_Y_noise_buffer_RCT,1)))); % idx is logical indices
            
            % 50,000 first line AED group
            generate_Y_noise_1st_Line = generate_Y_noise_buffer_RCT(Shuffle_indices(idx_sub_1st_Line),:);
            
            % 50,000 placebo group
            generate_Y_noise_buffer_Placebo_Control = generate_Y_noise_buffer_RCT(Shuffle_indices(~idx_sub_1st_Line),:);
            
            % Performing t test to make sure both first line AED and
            % Placebo has same IIC burden distribution between patients
            [h,p] = ttest2(sum(generate_Y_noise_buffer_Placebo_Control,2),sum(generate_Y_noise_1st_Line,2));
            
            
            
            
            if h==0
                
                
                
                % randomly generating attenuations for the 1 st line AED
                Trail_Drug_Affect  = rand(1,size(generate_Y_noise_1st_Line,1));
                %Trail_Drug_Affect  = 0.2*ones(1,size(generate_Y_noise_1st_Line,1));
                
                % randomly shuffling the population so some will have
                % complete cessation
                shuffle_population = randperm(size(generate_Y_noise_1st_Line,1));
                
                % Assume percentage of population has trail drug has 100%
                % effect throughout drug intervension
                
                Ind_Select_efficacy  = round(treatment_group_pop_sel*size(generate_Y_noise_1st_Line,1));
                Trail_Drug_Affect(shuffle_population(1:Ind_Select_efficacy))=0;
                % round decimal to two decimals
                % Trail drug effect is attenuation value
                Trail_Drug_Affect = round(Trail_Drug_Affect*100)/100;
                
                % Buffer to store outcomes of 50000 patients in first line
                % AED
                Final_Outcomes_1st_Line = zeros(size(generate_Y_noise_1st_Line,1),1);
                
                
                % RCT outcome for first line AED ARM/ Intervension ARM
                % Each loop is one patient and loop of 50,000 patients
                for i_1st_Line = 1:size(generate_Y_noise_1st_Line,1)
                    seizure_onset = find(generate_Y_noise_1st_Line(i_1st_Line,:)>0,1);
                    
                    Outcome_Measure = Drug_Outcomes(Seizure_burden_time_course,generate_Y_noise_1st_Line(i_1st_Line,seizure_onset:end),Trail_Drug_Affect(1,i_1st_Line),T_delay,T_drug_res,Out_Come,Monitoring_Offset,(No_min_Hr/Time_window));
                    Final_Outcomes_1st_Line(i_1st_Line) = Outcome_Measure(1);
                    
                end
                
                % Buffer to store outcomes of 50000 patients in control
                % group
                Final_Outcomes_control_Placebo = zeros(size(generate_Y_noise_buffer_Placebo_Control,1),1);
                
                
                % RCT out come for Placbo/ Control ARM
                % Each loop is one patients and loop of 50,000 patients
                for i_Placebo = 1:size(generate_Y_noise_buffer_Placebo_Control,1)
                    seizure_onset = find(generate_Y_noise_buffer_Placebo_Control(i_Placebo,:)>0,1);
                    
                    Outcome_Measure = Drug_Outcomes(Seizure_burden_time_course,generate_Y_noise_buffer_Placebo_Control(i_Placebo,seizure_onset:end),1-Drug_Effect_Placebo,T_delay,T_drug_res,Out_Come,Monitoring_Offset,(No_min_Hr/Time_window));
                    Final_Outcomes_control_Placebo(i_Placebo) = Outcome_Measure(1);
                    
                    
                    
                end
                
                
                % Calculation of the sample size
                % Mean and Variance of Trail AED
                Mean_Trail_AED = mean(Final_Outcomes_1st_Line);
                Var_Trail_AED = var(Final_Outcomes_1st_Line);
                
                
                % Mean and Variance of control group
                Mean_Trail_Control = mean(Final_Outcomes_control_Placebo);
                Var_Trail_Control = var(Final_Outcomes_control_Placebo);
                
                
                Sample_Size = (4*(Var_Trail_AED+Var_Trail_Control) *(1.96+0.842)^2)/(Mean_Trail_Control-Mean_Trail_AED)^2;
                
                Effect_Size = abs((Mean_Trail_Control-Mean_Trail_AED));
                
                Sample_Size_Buff( ite_n)=Sample_Size;
                Effective_Size_Buff( ite_n)=Effect_Size;
                
                %increment boost
                ite_n = ite_n+1;
            end
        end
        
        
        %%%%% Sample size %%%%%
        % 50 percentile
        Sample_size_Median(i_outcome,itr_time) = round(prctile(Sample_Size_Buff,50),2)
        Sample_size_Median_round(i_outcome,itr_time) = ceil(prctile(Sample_Size_Buff,50)/2)*2
        
        % 2.5 percentile
        Sample_size_Low(i_outcome,itr_time) = round(prctile(Sample_Size_Buff,2.5),2)
        Sample_size_Low_round(i_outcome,itr_time) = ceil(prctile(Sample_Size_Buff,2.5)/2)*2
        
        % 97.5 percentile
        Sample_size_Upper(i_outcome,itr_time) = round(prctile(Sample_Size_Buff,97.5),2)
        Sample_size_Upper_round(i_outcome,itr_time) = ceil(prctile(Sample_Size_Buff,97.5)/2)*2
        
        
        %%%%% Effective size %%%%%
        % 50 percentile
        Effect_Size_Median(i_outcome,itr_time) = round(prctile(Effective_Size_Buff,50),2)
        
        % 2.5 percentile
        Effect_Size_Low(i_outcome,itr_time) = round(prctile(Effective_Size_Buff,2.5),2)
        
        % 97.5 percentile
        Effect_Size_Upper(i_outcome,itr_time) = round(prctile(Effective_Size_Buff,97.5),2)
    end
end