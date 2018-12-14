clc;clear
close all

eV2m = 1239.84193e-9;
gamma = 9.6419e+03;
aw = 2.5;
xlamd = 0.03;

xlamds = (1+aw^2)/2/gamma^2*xlamd;
eph = eV2m/xlamds;

trange = 20e-15;
nslice = trange*3e8/xlamds