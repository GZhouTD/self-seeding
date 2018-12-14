clc;clear
close all
clight = 3e8;


p1 = [11.23e-15,4985/0.511];
p2 = [12.03e-15,4970/0.511];

d_gamma_d_t = (p2(2)-p1(2))/(p1(2)+p2(2))*2/(p2(1)-p1(1))

xlamd = 0.03;
xlamds = 1.1698e-09;
Kund_fel = 2.5;

a = (1+Kund_fel^2)/Kund_fel*xlamds/(clight*xlamd)*d_gamma_d_t