function result = process_sujet(dataSubject,optionList)
% PROCESS_SUJET process all manip for the subject data dataSubject
% result = process_sujet(dataSubject,optionList)
% IN
%      - dataSubject :  subject's structure
%      - optionList: array of pairs {'manipName',option, ...}
% OUT
%      - data       : struct with one field by processing

% find all loader functions in the directory

VERBOSE = 1;
procNameList = tools.processorList();

if nargin <2
    optionList = {};
end


% loop over manips
for iProc=1:numel(procNameList)
    procName = procNameList{iProc} ;
    % check if option
    idxPar = tools.find_string(procName, optionList);
    if isempty(idxPar)
        option = [];
    else
        option = optionList{idxPar+1};
        optionList(idxPar:idxPar+1) = [];
    end
    % load data
    try
        resultProc = process_sujet_model(dataSubject,procName,option) ;
    catch err
        resultProc = [];
        if VERBOSE, fprintf('%s (%s)\n',err.message, (procName)); end
    end
    result.(procName) = resultProc;
end
