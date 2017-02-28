
function [data] = load_data_choice(subjectDir, sessionList)
% LOAD_DATA_gripchoice extracts data for all choice manip of the MBB Battery
% [data] = load_data_choice(subjectDir, sessionList)
%
% IN
%       - subjectDir is the source directory
%       - sessionList is an optional argument to select session (by default, all sessions are extracted)
%           in "dir" order. ie. sessionList=[5 2] will extract data in the 5th and the 2d files listed in natural order (as session 1 and session2)
% OUT
%       - data containing one field by submanip. Each submanip contains:
%       -------------------------------------------------------------------------------------------------------------
%           * condition: experimental condition
%               - sessionNumber        : converted to 1 to nSession
%               - trialNumber          : in a given session
%               - leftItemNumber       : left item number (index, see rating manip to get actual items)
%               - rightItemNumber      : right item number (index, see rating manip to get actual items)
%               - leftRatingValue      : left item value (from rating)
%               - rightRatingValue     : right item value (from rating)
%               - differenceRatingValue: difference rating (right - left)
%       -------------------------------------------------------------------------------------------------------------
%           * behavior: depending on subject choie
%               - RT                   : reaction time
%               -  sideChoice          : choice side [1] right [-1] left
%       -------------------------------------------------------------------------------------------------------------
%           * misc
%               - session
%                   - fileName         : name of the loaded file


data=struct;

subManipList={'choiceR_s', 'choiceR_im_s', 'choiceP_s', 'choiceE_s'}; % to use as MANIP_NAME
submanipFieldName={'choiceReward', 'choiceRewardPicture', 'choicePunishment', 'choiceEffort'}; % to use as field name in data
nSubManip=length(subManipList);

if nargin<2, sessionList = []; end

for iSubmanip = 1:nSubManip
    
    %% LOAD FILES
    % ===========================================================================
    % result files should include manip 'manipName'
    MANIP_NAME = subManipList{iSubmanip};
    
    try
        [fileList, nSession] = tools.load_file(subjectDir, MANIP_NAME, sessionList);
        % concatenate sessions
        allData = vertcat(fileList.data);
        
        %% CONDITIONS
        % ===========================================================================
        subdata.condition.sessionNumber=[];
        for iSession = 1:nSession
            [nTrial,~]=size(fileList(iSession).data);
            subdata.condition.sessionNumber =   [subdata.condition.sessionNumber iSession * ones(1,nTrial)] ;
        end
        subdata.condition.trialNumber          = allData(:, 1)';
        subdata.condition.leftItemNumber       = allData(:, 3)';                         % left item number (index, see rating manip to get actual items)
        subdata.condition.rightItemNumber      = allData(:, 5)';                         % right item number (index, see rating manip to get actual items)
        subdata.condition.leftRatingValue      = allData(:, 4)';                         % left item value (from rating)
        subdata.condition.rightRatingValue     = allData(:, 6)';                         % right item value (from rating)
        subdata.condition.differenceRatingValue      = subdata.condition.rightRatingValue -  subdata.condition.leftRatingValue; % difference rating (right - left)
        
        %% BEHAVIOUR
        % ===========================================================================
        
        subdata.behavior.RT  = allData(:, 7)';                          % RT
        subdata.behavior.sideChoice  = allData(:, 8)';                  % choice side: [1] right [-1] left
        
        %% MISC
        % ===========================================================================
        for iSession = 1:nSession
            subdata.misc.session(iSession).fileName = fileList(iSession).fileName ;
        end
    catch
        subdata = [];
    end  
    data.(submanipFieldName{iSubmanip})= subdata ;  
end

end