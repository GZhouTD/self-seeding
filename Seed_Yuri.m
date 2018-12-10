%% -------------------------------------------------------------- %%
%                     Wake seed Yuri's theory                             %
%% -------------------------------------------------------------- %%
clear;
cd('E:\Phd thesis\Chuan paper\High efficiency FEL\detuning\self-seed\detuning\spectrum\1thstage')
load('SASEspec0.mat')
cd('E:\SLAC-2016\Self_seeding\Self-Seeding Simulation code\matlab\Seed generation Yuri\')
format long
global e_charge h_Plank c_speed
e_charge = 1.602176565e-19;     % charge unit[C]
h_Plank  = 6.62607004e-34;      % Plank constant [J-sec]
c_speed  = 299792458;           % speed of light[m/sec]

input_crystal_parameters;
crystal_struc.cry_thickness = cry_thickness;        % crystal thickness [m]
crystal_struc.cry_bragg = cry_bragg;                % bragg angle [deg]
crystal_struc.cry_asymmetry = cry_asymmetry;        % asymmetric angle[deg]
crystal_struc.pho_energy = pho_energy;              % photon energy [eV]
crystal_struc.ele_suscept0 = ele_suscept0;          % electric susceptibility
crystal_struc.ele_susceptH = ele_susceptH;          % electric susceptibility
crystal_struc.ele_susceptHbar = ele_susceptHbar;          % electric susceptibility



pho_energy = crystal_struc.pho_energy;
wavelength = h_Plank*c_speed/pho_energy/e_charge;                    % wavelength
w0 = 2*pi*c_speed/wavelength;
Omega = f;

%cd('E:\SLAC-2016\Self_seeding\Self-Seeding Simulation code\matlab\Response function of crystal')
[R001,R00,R0H] = Transmission(crystal_struc,Omega);
Ef= Esase.*R00;
lamda = flipud(c_speed./f);
energy0 = trapz(lamda,flipud(abs(Ef).^2));
Etimedomain = ifft(ifftshift(Ef));
amplitude = abs(Etimedomain); 
dtt = 1/(f(2)-f(1));
tt = dtt/length(f)*(0:length(f)-1);
energy1 = trapz(tt,amplitude.^2);
Ewake = Etimedomain*sqrt(energy0/energy1);
%    Atime=np.abs(Ewake)+Atime
figure;
% plot(f,abs(Ef).^2)
hold on;
plot(tt*1e15*0.3,abs(Ewake).^2,'r')
%hold on;
%plot((tt-tt(15600))*1e15*0.3,unwrap(angle(Ewake)),'r')
%plotyy((tt-tt(15600))*1e15*0.3,abs(Ewake).^2,(tt-tt(15600))*1e15*0.3,unwrap(angle(Ewake)))
cd('E:\SLAC-2016\Self_seeding\Self-Seeding\Two-stage self-seeding 4keV LCLS2\4keVtest\slot11-2stage-zsep10\Secondstage\Yuri approach\')
header1 = '? VERSION=2.0';
header2 = '? ZPOS          PRADO        ZRAYL       ZWAIST     PHASE';
zsep = 10;
xlamds=3.101155620672385e-10;
Eseed = Ewake(5911:14910);
Spower = abs(Eseed).^2;
Sphase = angle(Eseed);
z = (0:length(Spower)-1)'*zsep*xlamds;
zw = zeros(length(Spower),1);
zr = ones(length(Spower),1)*10;
out = [z Spower zr zw Sphase];
outputFile = 'Seed0.in';
fid = fopen(outputFile,'wt');
fprintf(fid,'%s\n',header1);
fprintf(fid,'%s\n',header2);
fprintf(fid,'%-18.11e        %-18.11e      %-18.11e       %-18.11e      %-18.11e\n',out');
fclose(fid);

