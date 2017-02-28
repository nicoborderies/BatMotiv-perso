function data = load_sujet_manip(manipName, subjectDir, sessionList)
% LOAD_SUJET_MANIP load data of the task manipName for the subject in subjectDir
% [ data ] = load_manip(manipName, dataDir, sessionList)
% IN
%      - manipName  : used to select the load_data_manipName
%      - subjectDir    : path to the parent directory containing the subjects' data
%      - sessionList: subset of sessions to load. Leave empty to load all sessions (default)
% OUT
%      - data       : vector of data structs s

% default arg
if nargin < 3, sessionList = []; end

% check if manip loader exists
loaderName = ['+loaders/load_data_' manipName];
if exist(loaderName,'file') ~= 2
        error('*** Cannot find the function %s, check the manipName.',loaderName);
end

% check if directory exists
if isdir(subjectDir)~=1
    error('*** Cannot find directory %s, check the path.',subjectDir)
end

% load datafile
% double-loop allowing to load data even if some fields are missing !
loader_fun = str2func(['loaders.load_data_' manipName]);
data=feval(loader_fun, subjectDir, sessionList);



