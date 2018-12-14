function [I,tbin] = current_plot_dpa(t,data,info,binsize)
cur = data.cur;
np = info.np;
dt = info.dt/3e8;
q = cur*dt/np;

queue = 1:length(t);
qq = q(ceil(queue/np));

tmin = min(t);
tmax = max(t);

tbin = linspace(tmin,tmax,binsize);
dtbin = mean(diff(tbin));
I = zeros(size(tbin));

for ii = 1:length(tbin)
    qi = qq((t>tbin(ii)-dtbin/2)&(t<=tbin(ii)+dtbin/2));
    I(ii) = sum(qi)./dtbin;
end
    
