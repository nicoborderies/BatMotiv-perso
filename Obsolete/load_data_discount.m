
function [data] = load_data_discount(subjectDir, sessionList)
% LOAD_DATA_discount extracts data for all choice manip of the MBB Battery
% [data] = load_data_discount(subjectDir, sessionList)
%
% IN
%       - subjectDir is the source directory
%       - sessionList is an optional argument to select session (by default, all sessions are extracted)
%           in "dir" order. ie. sessionList=[5 2] will extract data in the 5th and the 2d files listed in natural order (as session 1 and session2)
% OUT
%       - data containing one field by submanip. Each submanip contains:
%       -------------------------------------------------------------------------------------------------------------
%           * condition: experimental condition
%               - sessionNumber                            : converted to 1 to nSession
%               - trialNumber                              : in a given session
%               - immediateItemNumber                      : immediate item number (index, see in misc to see actual items)
%               - delayedItemNumber                        : delayed item number (index, see rating manip to get actual items)
%               - immediateItemRatingValue                 : rating of the immediate rating (from rating)
%               - delayedItemRatingValue                   : rating of the delayed item (from rating)
%               - delay                                    : delay (of the delayed item) in days
%               - delayedItemSide                          : side of the delayed iitem [1] right [-1] left
%               - differenceRatingValue                    : difference rating (delayed - immediate)
%       -------------------------------------------------------------------------------------------------------------
%           * behavior: depending on subject choie
%               - RT                   : reaction time
%               -  sideChoice          : choice side [1] right [-1] left
%               - isDelayedItemChoice  : choice type: [1] delayed item [0] immediate item
%       -------------------------------------------------------------------------------------------------------------
%           * misc
%               - session
%                   - fileName         : name of the loaded file
%                   - item List


data=struct;

subManipList={'discountR_s', 'discountR_im_s', 'discountP_s', 'discountE_s'}; % to use as MANIP_NAME
typeItemList={'rewardlist', 'rewardlist' 'punishlist', 'effortlist'};            % to get list of item
submanipFieldName={'discountReward', 'discountRewardPicture', 'discountPunishment', 'discountEffort'}; % to use as field name in data
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
        subdata.condition.trialNumber                = allData(:, 1)';
        subdata.condition.immediateItemNumber        = allData(:, 2)';                         % immediate item number (index, see rating manip to get actual items)
        subdata.condition.delayedItemNumber          = allData(:, 4)';                         % delayed item number (index, see rating manip to get actual items)
        subdata.condition.immediateItemRatingValue   = allData(:, 3)';                         % rating of the immediate rating
        subdata.condition.delayedItemRatingValue     = allData(:, 5)';                         % rating of the delayed item (from rating)
        subdata.condition.delay                      = allData(:, 6)';                         % delay (of the delayed item) in days
        subdata.condition.delayedItemSide            = 2*allData(:, 7)'-1;                     % side of the delayed iitem [1] right [-1] left
        subdata.condition.differenceRatingValue      = subdata.condition.delayedItemRatingValue -  subdata.condition.immediateItemRatingValue; % difference rating (delayed - immediate)
        
        %% BEHAVIOUR
        % ===========================================================================
        
        subdata.behavior.RT  = allData(:, 12)';                         % RT
        subdata.behavior.sideChoice  = 2*(allData(:, 10)'-114)-1;       % choice side: [1] right [-1] left
        subdata.behavior.isDelayedItemChoice  = allData(:, 11)';        % choice type: [1] delayed item [0] immediate item
        
        %% MISC
        % ===========================================================================
        for iSession = 1:nSession
            subdata.misc.session(iSession).fileName = fileList(iSession).fileName ;
        end
        if isfield(fileList(iSession), typeItemList{iSubmanip})
           subdata.misc.itemList =  getfield(fileList(iSession), typeItemList{iSubmanip});
        end

    catch
        subdata = [];
    end
    
    data.(submanipFieldName{iSubmanip})=subdata;
    
end
end