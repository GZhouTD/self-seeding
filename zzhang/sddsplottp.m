function [X,Y,Z,I,output] = sddsplottp(foldname,plotname,Q)
sddsfile = [foldname,'/',plotname,'.out'];
plainfile = [foldname,'/',plotname,'.pla'];
covert_command = ['sdds2plaindata ',sddsfile,' ',plainfile,' -col=t -col=p -sep=''  '''];
[~,output1] = system(covert_command);
fid = fopen(plainfile);
np = fscanf(fid,'%d',1);
pdata = fscanf(fid,'%f',np*2);
fclose(fid);
rm_command = ['rm ',plainfile];
[~,output2] = system(rm_command);
output = {output1,output2};
t = pdata(1:2:np*2);
p = pdata(2:2:np*2)*0.511;
t = t - mean(t);
[X,Y,Z,~,~,I] = contour_plot_current(t,p,300,300,0,Q);
