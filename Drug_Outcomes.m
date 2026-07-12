function [Outcome_Measure] = Drug_Outcomes(SB_time,SB,drug_effect,T_delay,delay_drug_res,Out_Come,Monitoring_Offset,Samples_in_1Hr)


% Treatment delay
t1 = (T_delay*Samples_in_1Hr)+1;

% delay in drug response to be 100% efficacy

drug_t2 = t1+(delay_drug_res*Samples_in_1Hr);
% time for drug effect to introduce

% Gradual increase drug effect from treat delay to drug delay
drug_effect_var = linspace(1,drug_effect,(delay_drug_res*Samples_in_1Hr)+2);
drug_effect_var = drug_effect_var(2:end-1);

% Effect of drug response from treatment delay to drug response delay
SB(t1:drug_t2-1) = drug_effect_var.*SB(t1:drug_t2-1);

% Time drug has complete effect
SB(drug_t2:end) = drug_effect*SB(drug_t2:end); % Complete seizure cessation


% Calculating total seizure burden
if Out_Come=="tSB"
    
    if rem(length(SB(1:end)),6)==0
        Outcome_Measure = sum(mean(reshape(SB(1:end),Samples_in_1Hr,length(SB(1:end))/Samples_in_1Hr)));
    else % adding the residual if total length of seizure burden is not divided by 6
        SB_residual = SB((end-rem(length(SB(1:end)),6))+1:end);
        SB_residual_cal = (mean(SB_residual)*(length(SB_residual)/Samples_in_1Hr));
        Outcome_Measure = sum(mean(reshape(SB(1:(end-rem(length(SB(1:end)),6))),Samples_in_1Hr,length(SB(1:(end-rem(length(SB(1:end)),6))))/Samples_in_1Hr)))+ SB_residual_cal;
    end
    
    
    
    
    % Calculating post intervention seizure burden for 72 Hrs
elseif Out_Come=="pSB72"
    Time_ = 72; % 24 hr
    %Post_Intervention = Samples_in_1Hr*Time; % i hr
    Post_Intervention = Samples_in_1Hr*(Time_);
    Outcome_Measure = sum(mean(reshape(SB(t1:t1+Post_Intervention-1),Samples_in_1Hr,length(SB(t1:t1+Post_Intervention-1))/Samples_in_1Hr)));
    
    % Calculating post intervention seizure burden for 24 Hrs
elseif Out_Come=="pSB24"
    Time_ = 24; % 24 hr
    %Post_Intervention = Samples_in_1Hr*Time; % i hr
    
    Post_Intervention = Samples_in_1Hr*(Time_);
    Outcome_Measure = sum(mean(reshape(SB(t1:t1+Post_Intervention-1),Samples_in_1Hr,length(SB(t1:t1+Post_Intervention-1))/Samples_in_1Hr)));
    
    
elseif Out_Come=="pPSB"
    
    
    % find the peak value
    [Peak_val,Peak_pos] = max(SB(t1:end));
    % find seizure burden for 24 hours from peak seizure burden
    Post_Intervention = Samples_in_1Hr*24;
    peak_start = t1+Peak_pos-1;
    Outcome_Measure = sum(mean(reshape(SB(peak_start:peak_start+Post_Intervention-1),Samples_in_1Hr,length(SB(peak_start:peak_start+Post_Intervention-1))/Samples_in_1Hr)));
    
end



end