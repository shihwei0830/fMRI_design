function [nTrials_est,mu] = get_nTrials_est(nT_run,t_range,tr,theta)
%
% get estimated number of trials given your desired run length, and jittering parameters
% e.g., [nTrials_est,mu] = get_nTrials_est(600,[2 12],2,4)
%
%

t_longFixTotal = 20; % long fixation at beginning and end
t_stim = 2;
t_feedback = 2;
expected_isi = theta;
expected_iti = theta;

%estimate trial number given run time, trial design, and jitter
expected_trialT = t_stim + t_feedback + expected_isi + expected_iti;
expected_nTrials = floor((nT_run-t_longFixTotal)/expected_trialT);

% Get mean of jitter given jitter range, mean of exponential, and tr
% mu will be different from theta because we are sampling from truncated
% exponential
[sample,histo,mu] = get_expo_sample(t_range,tr,theta,expected_nTrials);

% Use updated jitter mean to update expected isi and iti.
expected_isi = mu; 
expected_iti = mu;

% Update nTrials estimate
nTrials_est = floor((nT_run-t_longFixTotal)/(t_stim + t_feedback + expected_isi + expected_iti)); 

fprintf('\n Given run time %.2f sec, trial design, and average isi/iti %.2f\n',nT_run,theta);
fprintf('\n expected number of trials is %i\n',nTrials_est);


