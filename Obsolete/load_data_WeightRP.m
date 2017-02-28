function [data] = load_data_WeightRP(subjectDir, sessionList)
% LOAD_DATA_WEIGHTING extracts data for the Weight(RE,PE or RP) task of the MBB Battery
% [data] = load_data_weighting(subjectDir, sessionList)
%
% IN
%       - subjectDir is the source directory
%       - sessionList is an optional argument to select session (by default, all sessions are extracted)
%           in "dir" order. ie. sessionList=[5 2] will extract data in the 5th and the 2d files listed in natural order (as session 1 and session2)
%
%
% OUT
%       - data containing the following field
%       -------------------------------------------------------------------------------------------------------------
%           * condition: experimental condition
%               - sessionNumber        : converted to 1 to nSession
%               - trialNumber          : in a given session
% 
%               - benefitNumber        : from 1 to 24 (reference to the appropriate benefit list)
%               - benefitRating        : rating of benefit item, from 1 to
%                                        100 (reward for RE task, punsihment for PE task, reward for RP task)
% 
%               - costNumber           : from 1 to 24 (reference to the appropriate cost list)
%               - costRating           : rating of cost item, from 1 to
%                                        100 (effort for RE task, effort for PE task, punishment for RP task)
%       -------------------------------------------------------------------------------------------------------------
%           * behavior: depending on subject choice
%               - choice               : [0] no-go / [1] go
%               - choiceSide           : [-1] left press / [1] right press
%               - RT                   : reaction time
%       -------------------------------------------------------------------------------------------------------------
%           * misc
%               - session
%                   - type             : trade-off type (reward-effort / punishment-effort / reward-punishment)
%                   - rewardlist      : items used for the rewarding stimuli 
%                   - punishlist      : items used for the punishing stimuli 
%                   - effortlist      : items used for the effortfull stimuli  

%% LOAD FILES
% ===========================================================================
% result files should include manip 'manipName'
MANIP_NAME = 'WeightRP_s';

% load files in a struct
if nargin<2, sessionList = []; end
[fileList, nSession] = tools.load_file(subjectDir, MANIP_NAME, sessionList);

% concatenate sessions
allData = vertcat(fileList.data);

%% CONDITIONS
% ===========================================================================
data.condition.sessionNumber=[];
for iSession = 1:nSession
    [nTrial,~]=size(fileList(iSession).data);
    data.condition.sessionNumber =   [data.condition.sessionNumber iSession * ones(1,nTrial)] ;
end
data.condition.trialNumber = allData(:, 1)'  ;
data.condition.benefitNumber = allData(:, 2)' ;                              
data.condition.benefitRating    = allData(:, 3)'  ;                              
data.condition.costNumber = allData(:, 4)' ;              
data.condition.costRating =  allData(:, 5)' ;     


%% BEHAVIOUR
% ===========================================================================
data.behavior.choice  = allData(:, 7)';                                     
data.behavior.choiceSide = allData(:, 8)';                                   
data.behavior.RT      = allData(:, 6)';                                     % in sec

%% MISC
% ===========================================================================
% loop across sessions
for iSession = 1:nSession
    % filenames
    data.misc.session(iSession).fileName = fileList(iSession).fileName ;
    % get trade-off type
    if isfield(fileList(1, iSession),'rewardlist') && isfield(fileList(1, iSession),'effortlist')
         data.misc.type{iSession} = 'RE' ;
    elseif isfield(fileList(1, iSession),'punishlist') && isfield(fileList(1, iSession),'effortlist')
         data.misc.type{iSession} = 'PE' ;
    elseif isfield(fileList(1, iSession),'rewardlist') && isfield(fileList(1, iSession),'punishlist')
         data.misc.type{iSession} = 'RP' ;
    end

end
end



