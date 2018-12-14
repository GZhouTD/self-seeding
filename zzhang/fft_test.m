clc;clear
close all
Q = 189.4225e-12;
sddsfile = 'elegant/scan/BC1COL/col_3/s1/DL2END.out';
plainfile = 'elegant/scan/BC1COL/col_3/s1/DL2END.pla';
covert_command = ['sdds2plaindata ',sddsfile,' ',plainfile,' -col=t -col=p -sep=''  '''];
if ~exist(plainfile,'file')
    system(covert_command);
end
fid = fopen(plainfile);
np = fscanf(fid,'%d',1);
pdata = fscanf(fid,'%f',np*2);
fclose(fid);
t = pdata(1:2:np*2);
p = pdata(2:2:np*2);
t = t - mean(t);

lambda = 4;
zsep = 1;
Wz = load(['elegant/wake/wiggler_wake_6per_',num2str(lambda),'um_filter.txt']);
tbin = Wz(:,1);
Wf = Wz(:,2);
t1 = tbin(tbin>0);
Wf = Wf(tbin>-max(t1));
tbin = tbin(tbin>-max(t1));
Wf = Wf*(zsep/6)*1;
dt = tbin(2) - tbin(1);

N = hist(t,tbin);
I = Q/np*N/dt;

I_filter = filter1d(I,tbin,4e-15);
I_filter = abs(I_filter);

figure
plot(tbin/1e-15,I)
hold on
plot(tbin/1e-15,I_filter)