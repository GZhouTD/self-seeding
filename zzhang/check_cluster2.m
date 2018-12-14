function [sub_bullet,sub_oak] = check_cluster2(njob)

sub_bullet = 0;
sub_oak = 0;

[~,bjobs] = system('bjobs -q beamphysics');
strpos = strfind(bjobs,'zzhang');
if length(strpos) < njob+1
    sub_oak = 1;
    return;
end

[~,bjobs] = system('bjobs -q beamphysics-mpi');
strpos = strfind(bjobs,'zzhang');
if length(strpos) < njob+5
    sub_bullet = 1;
    return;
end
