clc;clear
% close all
lambda = 2;
Wz = load(['elegant/wake/wiggler_wake_6per_',num2str(lambda),'um.txt']);
tbin = Wz(:,1);
Wf = Wz(:,2);
t1 = tbin(tbin>0);
Wf = Wf(tbin>-max(t1));
tbin = tbin(tbin>-max(t1));


w1 = Wf(3845:3993);
w1 = smooth1d(w1,4);
Wf(3845:3993) = w1;

w1 = Wf(3788:3840);
w1 = smooth1d(w1,4);
Wf(3788:3840) = w1;


w1 = Wf(3691:3783);
w1 = smooth1d(w1,4);
Wf(3691:3783) = w1;


w1 = Wf(3625:3684);
w1 = smooth1d(w1,4);
Wf(3625:3684) = w1;


w1 = Wf(3533:3617);
w1 = smooth1d(w1,4);
Wf(3533:3617) = w1;

w1 = Wf(3463:3526);
w1 = smooth1d(w1,4);
Wf(3463:3526) = w1;

w1 = Wf(3376:3456);
w1 = smooth1d(w1,4);
Wf(3376:3456) = w1;

figure
% hold on
plot(tbin,Wf)

% filename = ['elegant/wake/wiggler_wake_6per_',num2str(lambda),'um_filter.txt'];
% fid = fopen(filename,'w');
% for j = 1:length(tbin)
%   fprintf(fid,'%12.8e  %12.8e\n',tbin(j),Wf(j));
% end
% fclose(fid);