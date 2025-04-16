function [sample,histo,mu] = get_expo_sample(t_range,tr,theta,n_sample)
%
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
cdf_pTrunc = cumsum(p_truncated);
%figure(1);clf
%plot(t(1:n_iti),p,'k-','linewidth',2);
%hold on;
%plot(t(1:n_iti),p_truncated,'b-','linewidth',2);
%plot(t(1:n_iti),cdf_pTrunc,'r','linewidth',2);
%legend original truncated cdf

%n_trials = 100;
%ef = n_trials*p_truncated
x=rand(n_sample,1);     % simulate from uniform [0,1]
histo = zeros(2,n_iti);
histo(1,:) = t(1:n_iti);
sample = zeros(n_sample,1);

for i=1:n_sample
    for j=1:n_iti
        if j==1
            if x(i)<=cdf_pTrunc(j)
                sample(i)=t(1);
                histo(2,j) = histo(2,j)+1;
            end
            
        else
            if x(i)<cdf_pTrunc(j) && x(i)>cdf_pTrunc(j-1)
                sample(i)=t(j);
                histo(2,j) = histo(2,j)+1;
            end 
           
        end
    end
end