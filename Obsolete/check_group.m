function s=check_group(group)

nSubject = numel(group.sujet);

s ='';

s=[s sprintf('\n  _____________________________________________________\n')];
s=[s sprintf('  = LOADING OPTIONS ===================================\n\n')];
s=[s sprintf('  +%19s : %s\n','root directory',group.misc.dataDir)];
s=[s sprintf('\n')];

s=[s sprintf('  +%19s : ','list of subjects')];
for i=1:numel(group.misc.subjectList)
    if i>1, s=[s sprintf('%25s','')]; end
    s=[s sprintf('%02d - %s \n',i,group.misc.subjectList{i})];
end
s=[s sprintf('\n')];

s=[s sprintf('  +%19s : ','excluded subjects')];
for i=0:numel(group.misc.excludedSubject)/2-1
    if i>1, s=[s sprintf('%25s','')]; end
    s=[s sprintf('%s - ',group.misc.excludedSubject{2*i+1})];
    s=[s sprintf('%02d ',group.misc.excludedSubject{2*i+2})];
    s=[s sprintf('\n')];
end
s=[s sprintf('\n')];

s=[s sprintf('  +%19s : ','excluded sessions')];
for i=0:numel(group.misc.excludedSession)/2-1
    if i>1, s=[s sprintf('%25s','')]; end
    s=[s sprintf('%s - ',group.misc.excludedSession{2*i+1})];
    s=[s sprintf('%02d ',group.misc.excludedSession{2*i+2})];
    s=[s sprintf('\n')];
end


s=[s sprintf('\n  _____________________________________________________\n')];
s=[s sprintf('  = TASKS SUMMARY =====================================\n\n')];

s=[s sprintf('%21s ','# sujet')];
s=[s sprintf(' %2d',1:nSubject)];
s=[s sprintf('\n')];
s=[s sprintf('      ---------------  ')];
s=[s sprintf(repmat('---',1,nSubject))];
s=[s sprintf('\n')];

manipList = fieldnames(group.sujet(1).data);
for iManip = 1:numel(manipList)
    manipName = manipList{iManip};
    s=[s sprintf('%21s ',manipName)];
    
    for iSujet = 1:nSubject
        if isstruct(group.sujet(iSujet).data.(manipName)) ;
            try
                nSession = num2str(numel(group.sujet(iSujet).data.(manipName).misc.session));
            catch
                nSession = '#';
            end
            s=[s sprintf('%3s',nSession)];
        else
            s=[s sprintf('  -')];
        end
    end
    s=[s sprintf('\n')];
end

if isfield(group.sujet(1),'result')
s=[s sprintf('\n  _____________________________________________________\n')];
s=[s sprintf('  = ANALYSES SUMMARY ==================================\n\n')];

s=[s sprintf('%21s ','# sujet')];
s=[s sprintf(' %2d',1:nSubject)];
s=[s sprintf('\n')];
s=[s sprintf('      ---------------  ')];
s=[s sprintf(repmat('---',1,nSubject))];
s=[s sprintf('\n')];

manipList = fieldnames(group.sujet(1).result);
for iManip = 1:numel(manipList)
    manipName = manipList{iManip};
    s=[s sprintf('%21s ',manipName)];
    
    for iSujet = 1:nSubject
        if isstruct(group.sujet(iSujet).result.(manipName)) ;
            s=[s sprintf('  x')];
        else
            s=[s sprintf('  .')];
        end
    end
    s=[s sprintf('\n')];
end
end

s=[s sprintf('\n  _____________________________________________________\n\n')];


end