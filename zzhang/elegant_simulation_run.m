function elegant_simulation_run(inp_struc)

sub_bullet = 0;
sub_oak = 0;

charge = inp_struc.charge;
dump_BC1BEG = inp_struc.dump_BC1BEG;
dump_BC1COL = inp_struc.dump_BC1COL;
dump_BC1MID = inp_struc.dump_BC1MID;
dump_BC1END = inp_struc.dump_BC1END;
dump_BC2END = inp_struc.dump_BC2END;
dump_L3END = inp_struc.dump_L3END;
dump_DL2END = inp_struc.dump_DL2END;
dump_XLBEG = inp_struc.dump_XLBEG;
dump_UNDBEG = inp_struc.dump_UNDBEG;
L1AMP = inp_struc.L1AMP;
L1PHASE = inp_struc.L1PHASE;
L2AMP = inp_struc.L2AMP;
L2PHASE = inp_struc.L2PHASE;
L1XPHASE = inp_struc.L1XPHASE;
L1XAMP = inp_struc.L1XAMP;
L3AMP = inp_struc.L3AMP;
L3PHASE = inp_struc.L3PHASE;
BC1ANG = inp_struc.BC1ANG;
BC2ANG = inp_struc.BC2ANG;
DL2R56 = inp_struc.DL2R56;
col_length = inp_struc.col_length;
x_max_col = inp_struc.x_max_col;
dx_col = inp_struc.dx_col;
parallelQ = inp_struc.parallelQ;


if parallelQ
    while 1
        [sub_bullet,sub_oak] = check_cluster2(1);
        if sub_bullet || sub_oak
            break;
        else
            pause(20);
        end
    end
end

foldname = inp_struc.main_fold;
temp_fold = [foldname,'/LCLS_scan_templete'];
temp_lte = [temp_fold,'/',inp_struc.lte_name];
temp_ele = [temp_fold,'/',inp_struc.ele_name];

% temp_lte = [temp_fold,'/LCLS25Aug2015W.lte'];
% temp_ele = [temp_fold,'/LCLS25Aug2015W.ele'];

scan_fold = [foldname,'/',inp_struc.scan_name];
if ~exist(scan_fold,'dir')
    mkdir(scan_fold);
end
[~,dirQ] = system(['ls ',scan_fold]);
ssnum = 1;
if isempty(dirQ)
    scan_run = [scan_fold,'/s1'];
    mkdir(scan_run);
else
    while 1
        fileQ = strfind(dirQ,['s',num2str(ssnum),' ']);
        if ~isempty(fileQ)
            ssnum = ssnum + 1;
        else
            while exist([scan_fold,'/s',num2str(ssnum)],'dir')
                ssnum = ssnum + 1;
            end
            scan_run = [scan_fold,'/s',num2str(ssnum)];
            mkdir(scan_run);
            break;
        end
    end
end

[A,Nline] = read_file(temp_lte);
for ii = 1:Nline
    A{ii} = strrep(A{ii},'$charge$',num2str(charge,'%e'));
    A{ii} = strrep(A{ii},'$dump_BC1BEG$',num2str(dump_BC1BEG));
    A{ii} = strrep(A{ii},'$dump_BC1COL$',num2str(dump_BC1COL));
    A{ii} = strrep(A{ii},'$dump_BC1MID$',num2str(dump_BC1MID));
    A{ii} = strrep(A{ii},'$dump_BC1END$',num2str(dump_BC1END));
    A{ii} = strrep(A{ii},'$dump_BC2END$',num2str(dump_BC2END));
    A{ii} = strrep(A{ii},'$dump_L3END$',num2str(dump_L3END));
    A{ii} = strrep(A{ii},'$dump_DL2END$',num2str(dump_DL2END));
    A{ii} = strrep(A{ii},'$dump_XLBEG$',num2str(dump_XLBEG));
    A{ii} = strrep(A{ii},'$dump_UNDBEG$',num2str(dump_UNDBEG));
    A{ii} = strrep(A{ii},'$L1AMP$',num2str(L1AMP));
    A{ii} = strrep(A{ii},'$L1PHASE$',num2str(L1PHASE));
    A{ii} = strrep(A{ii},'$L2AMP$',num2str(L2AMP));
    A{ii} = strrep(A{ii},'$L2PHASE$',num2str(L2PHASE));
    A{ii} = strrep(A{ii},'$L1XPHASE$',num2str(L1XPHASE));
    A{ii} = strrep(A{ii},'$L1XAMP$',num2str(L1XAMP));
    A{ii} = strrep(A{ii},'$L3AMP$',num2str(L3AMP));
    A{ii} = strrep(A{ii},'$L3PHASE$',num2str(L3PHASE));
    A{ii} = strrep(A{ii},'$BC1ANG$',num2str(BC1ANG));
    A{ii} = strrep(A{ii},'$BC2ANG$',num2str(BC2ANG));
    A{ii} = strrep(A{ii},'$DL2R56$',num2str(DL2R56,'%e'));
    A{ii} = strrep(A{ii},'$col_length$',num2str(col_length,'%e'));
    A{ii} = strrep(A{ii},'$x_max_col$',num2str(x_max_col,'%e'));
    A{ii} = strrep(A{ii},'$dx_col$',num2str(dx_col,'%e'));
end
lte_filename = [scan_run,'/',inp_struc.lte_name];
write_file(A,Nline,lte_filename);

[A,Nline] = read_file(temp_ele);
ele_filename = [scan_run,'/',inp_struc.ele_name];
write_file(A,Nline,ele_filename);


fidpara = fopen([scan_run,'/locallog'],'a+');
fprintf(fidpara,'charge:\t   %e\n',charge);
fprintf(fidpara,'dump_BC1END:\t   %d\n',dump_BC1END);
fprintf(fidpara,'dump_BC2END:\t   %d\n',dump_BC2END);
fprintf(fidpara,'dump_L3END:\t   %d\n',dump_L3END);
fprintf(fidpara,'dump_DL2END:\t   %d\n',dump_DL2END);
fprintf(fidpara,'dump_UNDBEG:\t   %d\n',dump_UNDBEG);
fprintf(fidpara,'L1AMP:\t   %f\n',L1AMP);
fprintf(fidpara,'L1PHASE:\t   %f\n',L1PHASE);
fprintf(fidpara,'L2AMP:\t   %f\n',L2AMP);
fprintf(fidpara,'L2PHASE:\t   %f\n',L2PHASE);
fprintf(fidpara,'L1XAMP:\t   %f\n',L1XAMP);
fprintf(fidpara,'L1XPHASE:\t   %f\n',L1XPHASE);
fprintf(fidpara,'L3AMP:\t   %f\n',L3AMP);
fprintf(fidpara,'L3PHASE:\t   %f\n',L3PHASE);
fprintf(fidpara,'BC1ANG:\t   %f\n',BC1ANG);
fprintf(fidpara,'BC2ANG:\t   %f\n',BC2ANG);
fprintf(fidpara,'DL2R56:\t   %e\n',DL2R56);
fprintf(fidpara,'col_length:\t   %f\n',col_length);
fprintf(fidpara,'x_max_col:\t   %f\n',x_max_col);
fprintf(fidpara,'dx_col:\t   %e\n',dx_col);
fclose(fidpara);
        
if parallelQ

timenow = clock;
if sub_bullet
    disp([num2str(timenow(4),'%d'),':',num2str(timenow(5),'%d'),':',num2str(round(timenow(6)),'%d'),': submit jobs to bullet...'])
    disp(['ssnum = ',num2str(ssnum)]);
    [~,b] = bsub_bullet_elegant(scan_run,32,inp_struc.ele_name);
    disp(b)
end
if sub_oak
    disp([num2str(timenow(4),'%d'),':',num2str(timenow(5),'%d'),':',num2str(round(timenow(6)),'%d'),': submit jobs to oak...'])
    disp(['ssnum = ',num2str(ssnum)]);
    [~,b] = bsub_oak_elegant(scan_run,32,inp_struc.ele_name);
    disp(b)
end
pause(20);

else
    cur_fold = pwd;
    cd(scan_run)
    bsub_command = ['elegant ',inp_struc.ele_name];
    [~,b] = system(bsub_command);
    cd(cur_fold);
    disp(b);
end