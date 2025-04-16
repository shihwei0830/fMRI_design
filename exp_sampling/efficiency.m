% efficiency
clear all

% An example that correlation between regressors hurt design efficiency
n_events=2;
n_t = 100;
mu = [2 4];
sigma_1 = [1 1.5;1.5 3];
sigma_2 = [1 0;0 3];
c = [1 -1;1 0;0 1];
w_vect = [1 1 1];

x_1 = mvnrnd(mu,sigma_1,n_t); % correlated regressors
eff_1 = get_efficiency(x_1,c,w_vect);

x_2 = mvnrnd(mu,sigma_2,n_t); % not correlated
eff_2 = get_efficiency(x_2,c,w_vect);
% efficiency without correlation is larger
fprintf('\nSimulate design matrix from multivariate normal distribution\n')
fprintf('\n Without correlation, efficiency is %.2f; with correlation, efficiency is %.2f\n',eff_2,eff_1); 

% Define time points of a run
% length of a run in sec
nT_run = 600;
orig_resl = 0.1; % original resolution in time (sec)
TR = 2;
t_grid = 0:orig_resl:nT_run;
nT_in_resl = length(t_grid);
t_run = zeros(1,nT_in_resl);

% Jitter parameters
t_range = [2 12];
tr = 2;
theta = 4;

% Get an estimate of the number of trials
nTrials_est = get_nTrials_est(nT_run,t_range,tr,theta);


% Simulate jittered ISI and ITI.
n_sim = 10000;
for i=1:n_sim
    t_isi(:,i) = get_expo_sample(t_range,tr,theta,n_trials);
    t_iti(:,i) = get_expo_sample(t_range,tr,theta,n_trials);
end

% Compute efficiency and select maximum efficiency
% regressor timeseries
% in every trial, stimulus is followed by feedback
tau = 2; delta=2;
params_hrf = [tau delta];
leave_best = 5;
best_eff = zeros(1,leave_best);
eff = zeros(1,n_sim);
count_best=0;
for i=1:n_sim
    % Create 3-column matrix for each regressor
    t_end = 0;
    for j=1:n_trials
        t_start = t_end;
        stim_onset = t_start;
        feedback_onset = stim_onset + t_stim + t_isi(j,i);
        t_stim_3Col(j,:) = [stim_onset t_stim 1];
        t_feedback_3Col(j,:) = [feedback_onset t_feedback 1];
        t_end = feedback_onset + t_feedback + t_iti(j,i);
    end
    
    % Generate neural timeseries
    neural_ts_stim = generate_neural_ts(t_stim_3Col,t_grid);
    neural_ts_feedback = generate_neural_ts(t_feedback_3Col,t_grid);
    neural_ts_stim = neural_ts_stim(1:length(t_grid));
    neural_ts_feedback = neural_ts_feedback(1:length(t_grid));
    
    % Generate BOLD timeseries
    [BOLD_ts_stim,t_grid_downSamp] = generate_BOLD_ts(neural_ts_stim,params_hrf,orig_resl,tr,t_grid);
    [BOLD_ts_feedback,t_grid_downSamp] = generate_BOLD_ts(neural_ts_feedback,params_hrf,orig_resl,tr,t_grid);
    
    % Compute efficiency
    x = [BOLD_ts_stim BOLD_ts_feedback]; %design efficiency
    c = [1 -1]; % contrast matrix
    w_vect = 1; % weight vector
    eff(i) = get_efficiency(x,c,w_vect);
    if i <= leave_best
        best_eff(i) = eff(i);
        best_t_stim_3Col(i,:,:) = t_stim_3Col;
        best_t_feedback_3Col(i,:,:) = t_feedback_3Col;
        best_t_isi(:,i) = t_isi(:,i);
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
            
            % Redefine best_feedback_3Col
            best_t_feedback_3Col_new(id_kept,:,:) = best_t_feedback_3Col(id_best_eff(id_kept),:,:);
            best_t_feedback_3Col_new(id_new,:,:) = t_feedback_3Col;
            best_t_feedback_3Col = best_t_feedback_3Col_new;

            % Redefine best isi
            best_t_isi_new(:,id_kept) = best_t_isi(:,id_best_eff(id_kept));
            best_t_isi_new(:,id_new) = t_isi(:,i);
            best_t_isi = best_t_isi_new;
            
            % Redefine best iti
            best_t_iti_new(:,id_kept) = best_t_iti(:,id_best_eff(id_kept));
            best_t_iti_new(:,id_new) = t_iti(:,i);
            best_t_iti = best_t_iti_new;
        end

    end
    
end
fname = 'exp3';
save(fname,'eff','best_eff','best_t_stim_3Col','best_t_feedback_3Col','best_t_isi','best_t_iti','id_best_eff_orig')

figure;
histogram(eff);
% figure;
% plot(t_grid_downSamp,BOLD_ts_stim,'k','linewidth',2);
% hold on;
% plot(t_grid_downSamp,BOLD_ts_feedback,'r','linewidth',2);
% 






