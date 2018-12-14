function [A,N_line]=read_file(str)
forigin = fopen(str);
nline=0;
while ~feof(forigin)
    tline=fgetl(forigin);
    nline=nline+1;
    A{nline}=tline;
end
fclose(forigin);
N_line=nline;