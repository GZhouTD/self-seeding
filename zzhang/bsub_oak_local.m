function [a,b] = bsub_oak_local(foldname,Ncore,filename)
if ~exist('filename','var')
    filename = 'genesis.in';
end
cur_fold = pwd;
dpos = strfind(filename,'.');
logfile = [filename(1:dpos),'log'];
cd(foldname);
bsub_command = ['bsub -a mympi -q beamphysics -n ',num2str(Ncore),' -o ', logfile,' genesis2_mpi ',filename];
% bsub_command = ['bsub -a mympi -q beamphysics -n ',num2str(Ncore),' -o ', logfile,' genesis_BP_FFP ',filename];
disp(bsub_command);
[a,b] = system(bsub_command);
cd(cur_fold);