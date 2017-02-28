function dataGroup = process_group(dataGroup, optionList)
% PROCESS_GROUP process all manip for all subjects data dataGroup
% result = process_group(dataGroup, optionList)
% IN
%      - dataGroup         : structure of group data
%      - optionList: array of pairs {'manipName',option, ...}
% OUT
%      - result    

if nargin < 2
    optionList = {};
end

nSubject = numel(dataGroup.sujet);

for iSubject=1:nSubject
    dataGroup.sujet(iSubject).result = process_sujet(dataGroup.sujet(iSubject),optionList);
end

disp(check_group(dataGroup));