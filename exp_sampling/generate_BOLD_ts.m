function [ts_downSamp,t_grid_downSamp] = generate_BOLD_ts(neural_ts,params,resel,tr,t_grid)
%
%
%
%

% HRF:
% Plot the HRF with these parameter values:
tau = params(1);
delta = params(2);
tmin=0; tmax=30;
t_HRF = tmin:resel:tmax; % in seconds
tshift = max(t_HRF-delta,0);
HIRF = (tshift/tau).^2 .* exp(-tshift/tau) / (2*tau);
%figure(3); clf;
% Plot it
%plot(t,HIRF);
%title('Hemodynamic Impulse Response Function')
%ylabel('Hemodynamic response')
%xlabel('Time (sec)')

% For every choice of parameters, the volume of the function is 1 (or very
% close to 1 - it would be exactly one if we sampled more finely than once
% per second):
%sum(HIRF)

%%% Computing the fMRI response from the neural activity and the
%%% hemodynamic impulse response

% Now we use convolution to transform the neuralActivity into an
% fmri signal:
ts = conv(neural_ts,HIRF);
ts = ts(1:length(neural_ts));
ts_downSamp = ts(1:tr/resel:end);
t_grid_downSamp = t_grid(1:tr/resel:end);

% figure;
% length(ts)
% length(t_grid)
% plot(t_grid,ts,'k','linewidth',2);
% 
% figure;
% plot(t_grid_downSamp,ts_downSamp,'k','linewidth',2);
% hold on;
% plot(t_grid_downSamp,ts_downSamp,'k.','markersize',20);

