function [X,Y,Z,I,output] = sddsplotxy(foldname,plotname,Q)
sddsfile = [foldname,'/',plotname,'.out'];
plainfile = [foldname,'/',plotname,'.pla'];
covert_command = ['sdds2plaindata ',sddsfile,' ',plainfile,' -col=x -col=y -sep=''  '''];
[~,output1] = system(covert_command);
fid = fopen(plainfile);
np = fscanf(fid,'%d',1);
pdata = fscanf(fid,'%f',np*2);
fclose(fid);
rm_command = ['rm ',plainfile];
[~,output2] = system(rm_command);
output = {output1,output2};
x = pdata(1:2:np*2);
y = pdata(2:2:np*2);
[X,Y,Z,~,~,I] = contour_plot_current(x,y,300,300,0,Q);