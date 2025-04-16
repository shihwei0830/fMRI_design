clear all

tr=2; 
tmin = 2;
tmax = 12;
t = tmin:tr:tmax+tr;
theta = 4;
n_iti = length(t)-1;
for i=1:n_iti
    p(i)=cdf('exp',t(i+1),theta) - cdf('exp',t(i),theta);
    %p(i)=pdf('exp',t(i),theta);
end
total_p = cdf('exp',tmax+tr,theta) - cdf('exp',tmin, theta);
p_truncated = p./total_p;
mu_t = sum(t(1:n_iti).*p_truncated)
cdf_pTrunc = cumsum(p_truncated);
figure(1);clf
plot(t(1:n_iti),p,'k-','linewidth',2);
hold on;
plot(t(1:n_iti),p_truncated,'b-','linewidth',2);
plot(t(1:n_iti),cdf_pTrunc,'r','linewidth',2);
legend original truncated cdf

n_trials = 100;
ef = n_trials*p_truncated
x=rand(n_trials,1);
nTrials_iti = zeros(1,n_iti);
for i=1:n_trials
    for j=1:n_iti
        if j==1
            if x(i)<=cdf_pTrunc(j)
                iti(i)=t(1);
                nTrials_iti(j) = nTrials_iti(j)+1;
            end
            
        else
            if x(i)<cdf_pTrunc(j) && x(i)>cdf_pTrunc(j-1)
                iti(i)=t(j);
                nTrials_iti(j) = nTrials_iti(j)+1;
            end 
           
        end
    end
end