function sujet = load_sujet(subjectDir, sessionListByManip)
% LOAD_SUJET load all manip for the subject data dataDir
% [ data ] = load_sujet(dataDir, sessionListByManip)
% IN
%      - dataDir    : path to the parent directory containing the subjects' data
%      - sessionListByManip: array of pairs {'manipName',sessionList, ...}
%                     indicating which sessions to load for each manip(default: all)
% OUT
%      - data       : struct with one field by manip + misc



% find all loader functions in the directory
manipNameList = tools.loaderList();

% check parameters
if nargin <2
    sessionListByManip = {};
end
for parName = sessionListByManip(1:2:end)
    if isempty(tools.find_string(parName{1}, manipNameList))
        error('*** %s is not a valid manip name, check the sessionListByManip parameters.',parName{1})
    end
end

% loop over manips
for iManip=1:numel(manipNameList)
    manipName = manipNameList{iManip} ;
    % check if session restriction
    idxPar = tools.find_string(manipName, sessionListByManip);
    if isempty(idxPar)
        sessionList = [];
    else
        sessionList = sessionListByManip{idxPar+1};
        sessionListByManip(idxPar:idxPar+1) = [];
    end
    % load data
    try
        dataManip = load_sujet_manip(manipName, subjectDir, sessionList) ;
    catch err
        dataManip = [];
        fprintf('%s\n',err.message);
    end
    data.(manipName) = dataManip;
end

sujet.data = tools.expandSubmanips(data);
sujet.misc.subjectDirectory =  subjectDir;


