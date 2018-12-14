clc;clear
close all

foldnum = 1;
Eph = 500;
gamma_modulation = 50;
clight = 3e8;
aw0 = 2.5;
xlambda = 2.58e-9;
lambdu = 3e-2;

[pdata,Q,npart] = read_distfile(['genesis/s',num2str(foldnum),'/beam_core_match.dist']);
z = pdata{5};
z = z - min(z);
gamma = pdata{6};

[zbin,gamma_avg,gamma_std] = sliceEnergy(z,gamma,300);

dz = diff(zbin);
z_avg = (zbin(2:end)+zbin(1:(end-1)))/2;
dg = diff(gamma_avg);
g_avg = (gamma_avg(2:end)+gamma_avg(1:(end-1)))/2;

dg_dt = dg./g_avg./(dz'/clight);

daw_dz_initial = (1+aw0^2)/aw0*xlambda/(clight*lambdu)*dg_dt;

zstart = 0;
zstep = 115;


figure;
yyaxis left
plot(z_avg/xlambda,daw_dz_initial);
hold on

for ii = 1:6
    plot([zstart+(ii-1)*zstep,zstart+(ii-1)*zstep],[-0.01,0.01],'b--')
end
hold off
xlabel('z (\lambda_s)')
ylabel('daw/dz')
enhance_plot('times',16,2,8)
legend off
yyaxis right
plot(z_avg/xlambda,g_avg);
xlabel('z (\lambda_s)')
ylabel('\gamma')
enhance_plot('times',16,2,8)
legend off