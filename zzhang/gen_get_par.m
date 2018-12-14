function pval = gen_get_par(infn,pstr,till_line_no)
%**pval = gen_get_par(infn,pstr)
%  To get parmeter from genesis input file or main output file;
%<<Input
%  infn 
%    - Genesis input/main output filename
%      Main output filename would be better.     
%  pstr 
%    - string for the name of parameter
%  till_line_no 
%    - stop searching the parameter at the till_line_no-th line
%      default: 165  
%>>Output
%  pval 
%    - value of the parameter
%**Modification log
%  eatablished on 2014-04-17-Thu (S.Huang)
%  modified on 2015-01-22-Thu (S.Huang)

if nargin<3||isempty(till_line_no)
  till_line_no = 165;
end
fid = fopen(infn,'r');
jj = 1;
while ~feof(fid) && jj<till_line_no
  tline = regexprep(fgetl(fid),' ','');
  if ~isempty(findstr(tline,[lower(pstr),'=']))||~isempty(findstr(tline,[upper(pstr),'=']))
    ideq = findstr(tline,'=');
    pval = eval(tline(ideq+1:end));
    break
  end
  jj = jj+1;
end
fclose(fid);

if ~exist('pval','var')
  disp('Warning: The parameter is not found ...');
  pval = [];
end
