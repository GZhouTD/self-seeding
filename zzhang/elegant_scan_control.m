clc;clear
close all

dump_BC1END = 1;
dump_BC2END = 0;
dump_L3END = 0;
dump_UNDBEG = 0;
L1AMP = 1;
L1PHASE = 0;
L2AMP = 1;
L2PHASE = 0;
L1XPHASE = 0;
L1XAMP = 1;
L3AMP = 1;
L3PHASE = 0;
BC1ANG = 1;
BC2ANG = 1;
DL2R56 = 0;

foldname = 'elegant/scan';
temp_fold = [foldname,'/LCLS_scan_templete'];
temp_lte = [temp_fold,'/LCLS25Aug2015W.lte'];
temp_ele = [temp_fold,'/LCLS25Aug2015W.ele'];

scan_fold = [foldname,'/LCLS_replace_test'];
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
    A{ii} = strrep(A{ii},'$dump_BC1END$',num2str(dump_BC1END));
    A{ii} = strrep(A{ii},'$dump_BC2END$',num2str(dump_BC2END));
    A{ii} = strrep(A{ii},'$dump_L3END$',num2str(dump_L3END));
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
    A{ii} = strrep(A{ii},'$DL2R56$',num2str(DL2R56));
end
lte_filename = [scan_run,'/LCLS25Aug2015W.lte'];
write_file(A,Nline,lte_filename);

[A,Nline] = read_file(temp_ele);
ele_filename = [scan_run,'/LCLS25Aug2015W.ele'];
write_file(A,Nline,ele_filename);

sub_bullet = 0;
sub_oak = 0;

while 1
    [sub_bullet,sub_oak] = check_cluster2(1);
    if sub_bullet || sub_oak
        break;
    else
        pause(20);
    end
end

fidpara = fopen([scan_run,'/locallog'],'a+');
fprintf(fidpara,'dump_BC1END:\t   %d\n',dump_BC1END);
fprintf(fidpara,'dump_BC2END:\t   %d\n',dump_BC2END);
fprintf(fidpara,'dump_L3END:\t   %d\n',dump_L3END);
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
fprintf(fidpara,'DL2R56:\t   %f\n',DL2R56);
fclose(fidpara);
        
        
timenow = clock;
if sub_bullet
    disp([num2str(timenow(4),'%d'),':',num2str(timenow(5),'%d'),':',num2str(round(timenow(6)),'%d'),': submit jobs to bullet...'])
    disp(['ssnum = ',num2str(ssnum)]);
    [~,b] = bsub_bullet_elegant(scan_run,128);
    disp(b)
end
if sub_oak
    disp([num2str(timenow(4),'%d'),':',num2str(timenow(5),'%d'),':',num2str(round(timenow(6)),'%d'),': submit jobs to oak...'])
    disp(['ssnum = ',num2str(ssnum)]);
    [~,b] = bsub_oak_elegant(scan_run,128);
    disp(b)
end