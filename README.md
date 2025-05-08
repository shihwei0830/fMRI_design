# fMRI_design
Use simulations to find optimal design efficiency for fMRI experiments

5/8/2025: For trial sequence where there is only one stimulus (no feedback followed) during a trial
get_optEff_1ev_deconv.m
get_optEff_1ev.m

example:
% n_trials = 49; nT_run = 520; t_range = [2 12]; tr = 2; mu = 5; n_sim = 10000; leave_best = 10; fname = 'exp0';
% results = get_optEff_1ev_deconv(n_trials,nT_run,t_range,tr,mu,n_sim,leave_best,fname)


4/16/2025: For trial sequence where a stimulus is always followed by feedback.
ISI is between stimulus and feedback; ITI is between feedback on trial n-1 and stimulus on trial n. Therefore, a task with sequential dependency.

Basic idea: sample ISI and ITI from truncated exponential distribution [minT maxT].  

To run simulations without constraining block/run time, run get_optimal_efficiency.
example: 
n_trials = 40;
t_range = [2 12; 2 12];  % first row for ISI; 2nd row for ITI; 1st column is minT, 2nd maxT.
tr = 2;
mu = [4 4];  % first element for ISI; 2nd for ITI
leave_best = 5;
fname = 'exp0';
results = get_optimal_efficiency(n_trials,t_range,tr,mu,n_sim,leave_best,fname);

To run simulations with constraining block/run time, run get_optimal_efficiency_conT.
example:
n_trials = 40;
t_range = [2 12; 2 12];  % first row for ISI; 2nd row for ITI
tr = 2;
mu = [4 4];  % first element for ISI; 2nd for ITI
leave_best = 5;
fname = 'exp0';
UB_yes = 1;
results = get_optimal_efficiency_conT(n_trials,t_range,tr,mu,n_sim,leave_best,fname,UB_yes);
