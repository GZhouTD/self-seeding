function [pdata,Q,npart] = read_distfile(filename)
fid = fopen(filename);
tstr1 = textscan(fid,'%s%s%s%s',2);
tstr2 = textscan(fid,'%s%s%s%f',2);
Q = tstr2{4}(1);
npart = tstr2{4}(2);
tstr3 = textscan(fid,'%s%s%s%s%s%s%s%s',1);
pdata = textscan(fid,'%f%f%f%f%f%f',npart);
fclose(fid);