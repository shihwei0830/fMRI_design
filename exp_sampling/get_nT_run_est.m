function nT_run_est = get_nT_run_est(mu,n_trials)
%
% get run length estimate given mean isi and its (mu) and number of trials
%
%

t_longFixTotal = 20; % long fixation at beginning and end
t_stim = 2;
t_feedback = 2;
expected_isi = mu;
expected_iti = mu;

nT_run_est = t_longFixTotal + ceil((t_stim + t_feedback + expected_isi + expected_iti)*n_trials);
