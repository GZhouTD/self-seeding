clc;clear
% close all

foldname = 'elegant/scan/BC1COL/col_3';
ssnum = 64;
zsep = 1;
nwigg = 18;
nwigg2 = 24;

x_max_col = 3;

charge = 189.4225e-12; % x_max_col = 3;

Q = charge;
plotname = 'XLBEG';
lambda = 2;
Wz = load(['elegant/wake/wiggler_wake_6per_',num2str(lambda),'um_filter.txt']);
tbin = Wz(:,1);
Wf = Wz(:,2);
t1 = tbin(tbin>0);
Wf = Wf(tbin>-max(t1));
tbin = tbin(tbin>-max(t1));
Wf = Wf*(zsep/6)*1;
dt = tbin(2) - tbin(1);

DL2R56 = 0e-6;
% R56_xleap1 = -390*1e-6;
% R56_xleap2 = -110*1e-6;
R56_xleap1 = -(200:10:500)*1e-6;
R56_xleap2 = -(0:10:300)*1e-6;

for ii = 1:length(ssnum)
    sddsfile = [foldname,'/s',num2str(ssnum),'/',plotname,'.out'];
    plainfile = [foldname,'/s',num2str(ssnum),'/',plotname,'.pla'];
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
    delta = (p-mean(p))/mean(p);
    t = t + DL2R56*delta/3e8;
    nsteps = nwigg/zsep;
    nsteps2 = nwigg2/zsep;
    R56 = -2*lambda*1e-6*zsep;
    T566 = -1.5*R56;
    
    for jj = 1:nsteps
        N = hist(t,tbin);
        I = Q/np*N/dt;
        I = filter1d(I,tbin,4e-15);
        I = abs(I);
        dE = conv(I*dt,Wf,'same');
        dE = dE';
        dp = interp1(tbin,dE/0.511,t);
        p = p + dp;
        delta = (p-mean(p))/mean(p);
        t = t + R56*delta/3e8+T566*delta.^2/3e8;
        t = t - mean(t);
    end
    t_bake = t;
    p_bake = p;
    output = cell(length(R56_xleap1),length(R56_xleap2));
    for s1 = 1:length(R56_xleap1)
        for s2 = 1:length(R56_xleap2)
            [s1,s2]
            t = t_bake;
            p = p_bake;
            delta = (p-mean(p))/mean(p);
            t = t + R56_xleap1(s1)*delta/3e8;
            [X1,Y1,Z1,~,~,I1] = contour_plot_current(t,p*0.511,300,300,0,Q);
            for jj = 1:nsteps2
                N = hist(t,tbin);
                I = Q/np*N/dt;
                I = filter1d(I,tbin,4e-15);
                I = abs(I);
                dE = conv(I*dt,Wf,'same');
                dE = dE';
                dp = interp1(tbin,dE/0.511,t);
                p = p + dp;
                delta = (p-mean(p))/mean(p);
                t = t + R56*delta/3e8+T566*delta.^2/3e8;
                t = t - mean(t);
                
            end
            [X2,Y2,Z2,~,~,I2] = contour_plot_current(t,p*0.511,300,300,0,Q);
            
            delta = (p-mean(p))/mean(p);
            t = t + R56_xleap2(s2)*delta/3e8;
            trange = [-40e-15,20e-15];
            prange = [4.92e3,5.01e3];
            [X3,Y3,Z3,~,~,I3] = contour_plot_current(t,p*0.511,600,600,0,Q);
%             output{s1,s2}.X1 = X1;
%             output{s1,s2}.Y1 = Y1;
%             output{s1,s2}.Z1 = Z1;
%             output{s1,s2}.I1 = I1;
%             output{s1,s2}.X2 = X2;
%             output{s1,s2}.Y2 = Y2;
%             output{s1,s2}.Z2 = Z2;
%             output{s1,s2}.I2 = I2;
            output{s1,s2}.X3 = X3;
            output{s1,s2}.Y3 = Y3;
            output{s1,s2}.Z3 = Z3;
            output{s1,s2}.I3 = I3;
        end
    end
end

figure
set(gcf,'position',[100,100,800,800])
subplot(2,1,1)
imagesc(X3*1e15,Y3/1000,Z3)
axis xy
colormap(jetvar)
xlabel('time (fs)')
ylabel('E (GeV)')
enhance_plot('times',16,2,8)
legend off
subplot(2,1,2)
plot(X3*1e15,I3)
xlabel('time (fs)')
ylabel('Current (kA)')
set(gca,'xlim',[X3(1)*1e15,X3(end)*1e15]);
set(gca,'ylim',[0, max(I3)*1.2]);
enhance_plot('times',16,2,8)
legend off