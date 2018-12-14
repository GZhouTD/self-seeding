function [a,b] = bsub_bullet_local(foldname,Ncore,filename)
if ~exist('filename','var')
    filename = 'genesis.in';
end
cur_fold = pwd;
dpos = strfind(filename,'.');
logfile = [filename(1:dpos),'log'];
cd(foldname);
% bsub_command = ['bsub -a mympi -q beamphysics-mpi -sla bpmpi -n ',num2str(Ncore),' -o ', logfile,' genesis_BP_FFP ',filename];
bsub_command = ['bsub -a mympi -q beamphysics-mpi -sla bpmpi -n ',num2str(Ncore),' -o ', logfile,' genesis2_mpi ',filename];
disp(bsub_command);
[a,b] = system(bsub_command);
cd(cur_fold);