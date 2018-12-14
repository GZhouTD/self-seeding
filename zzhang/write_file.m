function write_file(A,N_line,str)
    fout1=fopen(str,'wt');
    for nline=1:N_line
        fprintf(fout1,'%s',A{nline});
        fprintf(fout1,'\n');
    end
    fclose(fout1);