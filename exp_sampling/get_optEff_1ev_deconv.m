function results = get_optEff_1ev_deconv(n_trials,nT_run,t_range,tr,mu,n_sim,leave_best,fname)
%
% Single event (1eve) with deconvolution (FIR) 
%
% Trial sequence: cue (t_cue) --> choice (t_stim) --> feedback (t_feedback)
% --> ITI (t_iti)
%
%
% Total time (nT_run) as a constraint. It constrains the sum of ISI and ITI
% do not exceed some value, making sure that total run time is not
% exceeding nT_run.
%
% UB_yes: 1= add std of ITI distribution to the mean
% ITI for upper bound; 0 = mean of ITI distribution as
% upper bound
%
%
% example:
% n_trials = 48; nT_run = 540; t_range = [2 12]; tr = 2; mu = 5; n_sim = 20000; leave_best = 10; fname = 'exp0';
% results = get_optEff_1ev_deconv(n_trials,nT_run,t_range,tr,mu,n_sim,leave_best,fname)
%
%
% Run simulations to find optimal efficiency
% Before you run this,
% (1) Decide how long your run is in sec. e.g., 540 (sec); Say you want 48
% trials and it will take a total of 480 sec (including ITI). Set this
% higher than 480 because we have 10 sec fix before and 10 sec fix after,
% and we want to allow for more possible schedules to be included in
% selecting the best efficiency. In this example, set 540 sec to be total.
%
% (2) Set jittering parameters below
% e.g.
% t_range = [2 12]; % For setting up truncated exponential distribution on
% jitter time
% tr = 2;
% theta = 4; (mean of exponential distribution on jitter time)
%
% n_sim: number of simulations
% e.g., n_sim = 10000;
%
% leave_best: number of best efficiency to keep
% e.g., leave_best = 10;
%
% (3) fname: name of file you want to save your results. e.g., fname =
% 'exp0'
%
% Now you are ready to run this function.
%
%

% Define time points of a run
% length of a run in sec
orig_resl = 0.1; % original resolution in time (sec)

% Trial timing parameters
t_longFix_begin = 10;
t_longFix_end = 10;
t_longFixTotal = t_longFix_begin + t_longFix_end; % long fixation at beginning and end
t_stim = 4;
t_cue = 0.5;    % fixation cue to start the trial
t_feedback = t_cue; % feedback on the chosen option
t_event = t_cue + t_stim + t_feedback;
t_event_total = t_event*n_trials;
% Trial sequence: cue (t_cue) --> choice (t_stim) --> feedback (t_feedback)
% --> ITI (t_iti)

% Determine total ITI time
t_jitter_total = nT_run - t_longFixTotal - t_event_total;
t_iti_total = t_jitter_total;

% Simulate jittered ITI.
count_n_sim_iti = 0;
%t_extend = 60;
while 1
    t_iti_0 = get_expo_sample(t_range,tr,mu,n_trials);
    if sum(t_iti_0)<=t_iti_total %+ t_extend
        if count_n_sim_iti>0
            match_check = t_iti == t_iti_0;
            sum_match_check = sum(match_check);
            if ~any(sum_match_check == n_trials)
                count_n_sim_iti = count_n_sim_iti+1;
                t_iti(:,count_n_sim_iti) = t_iti_0;
                if count_n_sim_iti == n_sim
                    break
                end
            end
        else
            count_n_sim_iti = count_n_sim_iti+1;
            t_iti(:,count_n_sim_iti) = t_iti_0;
        end
    end
end


% Compute efficiency and select maximum efficiency
% regressor timeseries
% in every trial, stimulus is followed by feedback
tau = 2; delta=2;
params_hrf = [tau delta];
best_eff = zeros(1,leave_best);
eff = zeros(1,n_sim);
max_t_forDeconv=20;
t_forDeconv = 0:tr:max_t_forDeconv;
nT_deconv = length(t_forDeconv);
for i=1:n_sim
    % Create 3-column matrix for each regressor
    %t_end = t_longFix; % there is a 10s fixation before 1st trial
    t_end = t_longFix_begin;
    for j=1:n_trials
        t_start = t_end;
        stim_onset = t_start + t_cue;
        t_stim_3Col(j,:) = [stim_onset t_stim 1];
        t_end = t_start + t_event + t_iti(j,i);
    end
    
    % Generate neural timeseries
    nT_simRun = ceil(sum(t_iti(:,i))) + t_event_total + t_longFixTotal;
    t_grid = 0:orig_resl:nT_simRun;
    BOLD_ts_stim=[];
    % Set deconvolve (FIR) design matrix
    for k=1:nT_deconv
        t_stim_3Col_0 = t_stim_3Col;
        t_stim_3Col_0(:,1) = t_stim_3Col(:,1) + t_forDeconv(k);
        t_stim_3Col_0(:,2) = tr;    % duration of each FIR
        neural_ts_stim = generate_neural_ts(t_stim_3Col_0,t_grid);
        neural_ts_stim = neural_ts_stim(1:length(t_grid));
        
        % Generate BOLD timeseries
        [BOLD_ts_stim_0,t_grid_downSamp] = generate_BOLD_ts(neural_ts_stim,params_hrf,orig_resl,tr,t_grid);
        %size(BOLD_ts_stim_0);
        BOLD_ts_stim(:,k) = BOLD_ts_stim_0;
    end
    
    % Compute design efficiency
    x = BOLD_ts_stim; %design matrix
    c = eye(nT_deconv); % contrast matrix
    w_vect = ones(1,nT_deconv); % weight vector
    eff(i) = get_efficiency(x,c,w_vect);
    if i <= leave_best
        best_eff(i) = eff(i);
        best_t_stim_3Col(i,:,:) = t_stim_3Col;
        best_t_iti(:,i) = t_iti(:,i);
        if i == leave_best
            [best_eff,id_best_eff_orig] = sort(best_eff,'descend');
            % update best_eff according to sorted best_eff
            % orig: id in original simulation order, e.g., [2 1 3 4 5] indicates that
            % the biggest eff is from the second simulation, the second
            % biggest from the first simulation, etc.
        end
    else
        
        % Create test bed
        best_eff_test = [best_eff eff(i)];  % a new eff vector containing the best eff so far and the new eff;
        id_best_eff_orig_test = [id_best_eff_orig i];
        % e.g., say i=26, then id_best_eff_test = [2 1 3 4 5 26]
        
        % Sort the new eff vector to get the index from best to worst
        [best_eff_test_sorted,id_best_eff_test_sorted] = sort(best_eff_test,'descend');
        % id_best_eff_test_sorted will be numbers from 1 to 6, e.g., [3 4 1
        % 2 6 5] indicates that the 3rd element in best_eff_test has the
        % biggest eff, then the 4th element, and etc.
        
        % Update best eff
        best_eff = best_eff_test_sorted(1:leave_best);      % get the best eff from the sorted new eff vector
        
        % Update id_best_eff (always 5 of the numbers from 1 to 6)
        id_best_eff = id_best_eff_test_sorted(1:leave_best);
        
        % Update id_best_eff_orig
        id_best_eff_orig = id_best_eff_orig_test(id_best_eff_test_sorted(1:leave_best));    % get the id of the new best eff
        % e.g., id_best_eff = [1 26 6 4 16];
        
        
        % say xxx
        % if one of the best eff is the new eff,
        if any(id_best_eff == (leave_best+1))
            id_kept = find(id_best_eff ~= (leave_best+1));   % id_kept = [1 3 4 5]
            id_new = find(id_best_eff == (leave_best+1));    % e.g., the new eff is ranked 2nd, i.e., id_new = 2
            
            % Redefine best_stim_3Col
            best_t_stim_3Col_new(id_kept,:,:) = best_t_stim_3Col(id_best_eff(id_kept),:,:);
            % id_kept = [1 3 4 5];
            % id_best_eff(id_kept) = [1 6 4 16]
            best_t_stim_3Col_new(id_new,:,:) = t_stim_3Col;
            best_t_stim_3Col = best_t_stim_3Col_new;
            
            % Redefine best iti
            best_t_iti_new(:,id_kept) = best_t_iti(:,id_best_eff(id_kept));
            best_t_iti_new(:,id_new) = t_iti(:,i);
            best_t_iti = best_t_iti_new;
        end
        
    end
    
end

sum_iti = ceil(sum(best_t_iti));
sum_run_totalT = t_longFixTotal + t_event_total + sum_iti;
results.eff = eff;
results.best_eff = best_eff;
results.best_t_stim_3Col = best_t_stim_3Col;
results.best_t_iti = best_t_iti;
results.id_best_eff_orig = id_best_eff_orig;
results.sum_run_totalT = sum_run_totalT;
results.t_range = t_range;
results.mu = mu;
results.tr = tr;
results.n_trials = n_trials;
results.t_cue = t_cue;
results.t_stim = t_stim;
results.t_feedback = t_feedback;
results.t_longFixTotal = t_longFixTotal;

save(fname,'results');

%fprintf('\n run length from the best efficiency is %d\n',sum_run_totalT)

figure;
subplot(2,1,1)
histogram(eff);
title('efficiency');

subplot(2,1,2);
histogram(best_t_iti);
title('best ITI');

% figure;
% plot(t_grid_downSamp,BOLD_ts_stim,'k','linewidth',2);
% hold on;
% plot(t_grid_downSamp,BOLD_ts_feedback,'r','linewidth',2);
%






