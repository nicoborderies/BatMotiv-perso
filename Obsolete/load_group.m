function group = load_group(dataDir, subjectExcludedListByManip, sessionListByManip)
% LOAD_GROUP load all manip for all subjects in dataDir
% [ data ] = load_group(dataDir, subjectListByManip, sessionListByManip)
% IN
%      - dataDir           : path to the parent directory containing the subjects folders
%      - subjectExcludedListByManip  
%                          : array of {manipName,subjectList,...} pairs (default = all subjects)
%                            indicating which subjects to exclude in each manip
%      - sessionListByManip: array of {manipName,sessionList,...} pairs (default = all sessions)
% OUT
%      - data              : vector of data containing a 1*nSubject field [manipName] by manip


if nargin < 2
    subjectExcludedListByManip = {};
end
if nargin < 3
    sessionListByManip = {};
end

% find all subjects
subjectDirList = dir([dataDir filesep 'sub*']) ; 
subjectDirList = tools.sort_nat({subjectDirList.name}) ;

% find subject folders
if isempty(subjectDirList)
    error('*** Cannot find any subject folders in directory %s, check the path.',dataDir);
end

nSubjects = numel(subjectDirList);

% load data
for iSubject=1:nSubjects
    subjectDir = [dataDir filesep subjectDirList{iSubject}];
    group.sujet(iSubject) = load_sujet(subjectDir, sessionListByManip);
    group.sujet(iSubject).misc.number = iSubject;
end

% exclusion
manipNameList = fieldnames(group.sujet(1).data) ;
for iManip=1:numel(manipNameList)
    manipName = manipNameList{iManip} ;
    % check if subject restriction
    idxPar = tools.find_string(manipName, subjectExcludedListByManip);
    if ~isempty(idxPar)
        excludedSubjects = subjectExcludedListByManip{idxPar+1};
        if nSubjects < max(excludedSubjects)
           error('*** The directory %s only contains %d subjects, but you asked to exclude subjects: %s from the manip ''%s''.',subjectDir,nSubjects,num2str(excludedSubjects),manipName);
        end
        for iSubject = excludedSubjects
           group.sujet(iSubject).data.(manipName) = [];
        end
    end    
end

%%
group.misc.dataDir = dataDir ;
group.misc.subjectList = subjectDirList ;
group.misc.excludedSubject = subjectExcludedListByManip ;
group.misc.excludedSession = sessionListByManip ;

disp(check_group(group));



