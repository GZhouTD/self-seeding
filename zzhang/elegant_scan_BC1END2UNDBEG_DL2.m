clc;clear
close all
dump_BC1BEG = 1;
dump_BC1COL = 1;
dump_BC1MID = 1;
dump_BC1END = 1;
dump_BC2END = 0;
dump_L3END = 0;
dump_DL2END = 0;
dump_XLBEG = 0;
dump_UNDBEG = 0;
L1AMP = 1;
L1PHASE = 0;
L2AMP = 1;
L2PHASE = [0.25,0.5];
L1XPHASE = 0;
L1XAMP = 1;
L3AMP = 4.2;     % L3AMP=1 : 3.5GeV   L3AMP=4.2 : 5GeV
L3PHASE = 0;
BC1ANG = 1;
BC2ANG = 1;
DL2R56 = 200*1e-6;
col_length = 0.05;
dx_col = -0;
parallelQ = 1;

beam_num = 3;
x_max_col = beam_num*1e-3;

% charge = 249.2432e-12;  % x_max_col = 5;
% charge = 231.9349e-12; % x_max_col = 4;
% charge =  217.5695e-12; % x_max_col = 3.6;
charge = 189.4225e-12; % x_max_col = 3;


lte_name = 'BEAM2UNDBEG.lte';
ele_name = 'BEAM2UNDBEG.ele';

% lte_name = 'LCLS25Aug2015W.lte';
% ele_name = 'LCLS25Aug2015W.ele';

% lte_name = 'BEAM2BC1BEG.lte';
% ele_name = 'BEAM2BC1BEG.ele';

for jj = 1:length(L2PHASE)
    for ii = 1:length(DL2R56)
        inp_struc.charge = charge;
        inp_struc.dump_BC1BEG = dump_BC1BEG;
        inp_struc.dump_BC1COL = dump_BC1COL;
        inp_struc.dump_BC1MID = dump_BC1MID;
        inp_struc.dump_BC1END = dump_BC1END;
        inp_struc.dump_BC2END = dump_BC2END;
        inp_struc.dump_L3END = dump_L3END;
        inp_struc.dump_DL2END = dump_DL2END;
        inp_struc.dump_XLBEG = dump_XLBEG;
        inp_struc.dump_UNDBEG = dump_UNDBEG;
        inp_struc.L1AMP = L1AMP;
        inp_struc.L1PHASE = L1PHASE;
        inp_struc.L2AMP = L2AMP;
        inp_struc.L2PHASE = L2PHASE(jj);    
        inp_struc.L1XPHASE = L1XPHASE;
        inp_struc.L1XAMP = L1XAMP;
        inp_struc.L3AMP = L3AMP;
        inp_struc.L3PHASE = L3PHASE;
        inp_struc.BC1ANG = BC1ANG;
        inp_struc.BC2ANG = BC2ANG;
        inp_struc.DL2R56 = DL2R56(ii);
        inp_struc.col_length = col_length;
        inp_struc.x_max_col = x_max_col;
        inp_struc.dx_col = dx_col;
        inp_struc.lte_name = lte_name;
        inp_struc.ele_name = ele_name;
        inp_struc.parallelQ = parallelQ;
        inp_struc.main_fold = 'elegant/scan';
        inp_struc.scan_name = ['BC1COL/col_',num2str(beam_num)];
        
        elegant_simulation_run(inp_struc);
        
    end
end