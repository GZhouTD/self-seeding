%% -------------------------------------------------------------- %%
%                 The ratio of SASE power to seed power                   %
%% -------------------------------------------------------------- %%
function Spower = Seedpower(Ppulse,Undulength,LG3D,Rho,Lpulse,d,Ep,theta,chi0)
% ----------------------------------------------------------------------- %
% 
% INPUT:
%  Undulength, The Undulator length before the crystal  [m]
%  LG3D,  3-D gain length [m]
%  Rho,   Pierce parameter
%  Lpulse,Electron beam length
%  d,     The thickness of crystal, [m]
%  Ep,    photon energy, [eV]
%  theta, Bragg ange, [deg]
%  extin, Extinction length, [m]
%  chi0,  Electric susceptibility of the material
% ----------------------------------------------------------------------- %
% d=30e-6;
% Ep=4000;
% theta=48.924;
% extin=2.8688267611461509e-6;
% chi0 = 0.33863e-04 +1i*0.53162e-06;  %diamond electric susceptibility (111)
% chi0 = 0.66596e-5+1i*0.19046e-6;    %silicon electric susceptibility (111)

h = 6.62607004e-34;
a = 1.602176565e-19;
c = 299792458; 
l = h*c/Ep/a;           % wavelength
w0 = 2*pi*c/l;
K0 = w0/c;
NG = Undulength/LG3D;
extin = 2*pi*sind(theta)/(K0*abs(chi0));
Spower = Ppulse*0.0058*sqrt(NG/(2*pi))*l*Lpulse/(6*Rho)*(pi^2*d/(extin^2*sind(theta)))^2;
