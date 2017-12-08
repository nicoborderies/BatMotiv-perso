%% make_motiscan_mph_sess0_report

% data folder
datadir = 'B:\nicolas.borderies\projets\batmotiv\données\MPH';
cd(datadir);

% subject list
flist = cellstr(ls);
parser = @(str) contains(str,'sub');
select = cell2mat(cellfun(parser,flist,'UniformOutput',0));
sublist = flist(select);
 
% parameters
sess = 0; % training session
workingdir = 'B:\nicolas.borderies\projets\batmotiv\motiscan-battery\scripts';
cd(workingdir);

% subject loop
for isub = 1:numel(sublist)
   subid = str2num(sublist{isub}(4:end));
   check_behavior(subid,sess,datadir,'all')
end


