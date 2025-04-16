function [sample,histo,mu] = get_expo_sample_representative(t_range,tr,theta,n_sample)
%
% !!!Representative sampling
% Sample from truncated exponential
% t_range = [tmin tmax]
% tr (TR)
% theta: mean of exponential
% n_sample: sample size
%
% Outputs:
% sample: sample of the truncated exponential
% histo: frequency distribution of the sample
% mu: mean of the truncated exponential
% 
%

% Set up random seed for random number generator
rseed = sum(1000*clock);
rng(rseed)

% Compute truncated exponential
tmin = t_range(1);
tmax = t_range(2);
t = tmin:tr:tmax+tr;
n_iti = length(t)-1;
for i=1:n_iti
    p(i)=cdf('exp',t(i+1),theta) - cdf('exp',t(i),theta);
end
total_p = cdf('exp',tmax+tr,theta) - cdf('exp',tmin, theta);
p_truncated = p./total_p;
mu = sum(t(1:n_iti).*p_truncated);

histo(1,:) = t(1:n_iti);
histo(2,:) = ceil(p_truncated.*n_sample);

sample=[];
for i=1:n_iti
    sample = [sample; histo(1,i)*ones(histo(2,i),1)];
end

id = randsample(n_sample,n_sample);
sample = sample(id);

