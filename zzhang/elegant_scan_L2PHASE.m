clc;clear
close all

dump_BC1END = 1;
dump_BC2END = 0;
dump_L3END = 0;
dump_DL2END = 0;
dump_UNDBEG = 0;
L1AMP = 1;
L1PHASE = 0;
% L2AMP = 1;
L2PHASE = -(0.25:0.25:1.5);
L1XPHASE = 0;
L1XAMP = 1;
L3AMP = 1;
L3PHASE = 0;
BC1ANG = 1;
BC2ANG = 1;
DL2R56 = 0;

p_BC1END = 4.305163195e2;
p_BC2BEG = 5.87558e3;
e_L2_MeV = (p_BC2BEG - p_BC1END)*0.511;
nominal_phase_2 = -37.5;
amp_full_0 = e_L2_MeV/cosd(nominal_phase_2);
amp_full = e_L2_MeV./cosd(nominal_phase_2+L2PHASE);
L2AMP = amp_full./amp_full_0;


for ii = 1:length(L2PHASE)
    
    inp_struc.dump_BC1END = dump_BC1END;
    inp_struc.dump_BC2END = dump_BC2END;
    inp_struc.dump_L3END = dump_L3END;
    inp_struc.dump_DL2END = dump_DL2END;
    inp_struc.dump_UNDBEG = dump_UNDBEG;
    inp_struc.L1AMP = L1AMP;
    inp_struc.L1PHASE = L1PHASE;
    inp_struc.L2AMP = L2AMP(ii);
    inp_struc.L2PHASE = L2PHASE(ii);   
    inp_struc.L1XPHASE = L1XPHASE;
    inp_struc.L1XAMP = L1XAMP;
    inp_struc.L3AMP = L3AMP;
    inp_struc.L3PHASE = L3PHASE;
    inp_struc.BC1ANG = BC1ANG;
    inp_struc.BC2ANG = BC2ANG;
    inp_struc.DL2R56 = DL2R56;
    inp_struc.main_fold = 'elegant/scan';
    inp_struc.scan_name = 'L2PHASE';
    
    elegant_simulation_run(inp_struc);
    
end