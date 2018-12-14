clc;clear
close all
dump_BC1BEG = 0;
dump_BC1COL = 0;
dump_BC1MID = 0;
dump_BC1END = 0;
dump_BC2END = 0;
dump_L3END = 0;
dump_DL2END = 0;
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
DL2R56 = 0*1e-6;
col_length = 0.05;
x_max_col = 3.0e-3;
dx_col = -0;
parallelQ = 0;

lte_name = 'BC1.lte';
ele_name = 'BC1.ele';

% lte_name = 'LCLS25Aug2015W.lte';
% ele_name = 'LCLS25Aug2015W.ele';

% lte_name = 'BEAM2BC1BEG.lte';
% ele_name = 'BEAM2BC1BEG.ele';

for ii = 1:length(x_max_col)
    inp_struc.dump_BC1BEG = dump_BC1BEG;
    inp_struc.dump_BC1COL = dump_BC1COL;
    inp_struc.dump_BC1MID = dump_BC1MID;
    inp_struc.dump_BC1END = dump_BC1END;
    inp_struc.dump_BC2END = dump_BC2END;
    inp_struc.dump_L3END = dump_L3END;
    inp_struc.dump_DL2END = dump_DL2END;
    inp_struc.dump_UNDBEG = dump_UNDBEG;
    inp_struc.L1AMP = L1AMP;
    inp_struc.L1PHASE = L1PHASE;
    inp_struc.L2AMP = L2AMP;
    inp_struc.L2PHASE = L2PHASE;    inp_struc.L1XPHASE = L1XPHASE;
    inp_struc.L1XAMP = L1XAMP;
    inp_struc.L3AMP = L3AMP;
    inp_struc.L3PHASE = L3PHASE;
    inp_struc.BC1ANG = BC1ANG;
    inp_struc.BC2ANG = BC2ANG;
    inp_struc.DL2R56 = DL2R56;
    inp_struc.col_length = col_length;
    inp_struc.x_max_col = x_max_col(ii);
    inp_struc.dx_col = dx_col;
    inp_struc.lte_name = lte_name;
    inp_struc.ele_name = ele_name;
    inp_struc.parallelQ = parallelQ;
    inp_struc.main_fold = 'elegant/scan';
    inp_struc.scan_name = 'BC1MID';
    
    elegant_simulation_run(inp_struc);
    
end