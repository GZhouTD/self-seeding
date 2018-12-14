function [a,b] = bsub_oak_elegant(foldname,Ncore,filename)
if ~exist('filename','var')
    filename = 'LCLS25Aug2015W.ele';
end
cur_fold = pwd;
logfile = 'run.log';
cd(foldname);
bsub_command = ['bsub -a mympi -q beamphysics -n ',num2str(Ncore),' -o ', logfile,' Pelegant ',filename];
disp(bsub_command);
[a,b] = system(bsub_command);
cd(cur_fold);