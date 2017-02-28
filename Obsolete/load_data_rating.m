
function [data] = load_data_rating(subjectDir, sessionList)
% LOAD_DATA_rating extracts data for all rating manip of the MBB Battery
% [data] = load_data_rating(subjectDir, sessionList)
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
%               - itemNumber           : Item Number
%               - itemSubtype             : Item type: [0/1]
%                                        * Reward :      [0] = food; [1] = other
%                                        * Punishement : ???
%                                        * Effort :      [0] = physical; [1] = other
%       -------------------------------------------------------------------------------------------------------------
%           * behavior: depending on subject choice
%               - rating               : rating (betwenn 0 and 100)
%       -------------------------------------------------------------------------------------------------------------
%           * misc
%               - session
%                   - fileName         : name of the loaded file
%               - itemList             : item list


data=struct;

subManipList={'ratingR_s', 'ratingR_im_s', 'ratingP_s', 'ratingE_s'}; % to use as MANIP_NAME
typeItemList={'rewardlist', 'rewardlist' 'punishlist', 'effortlist'};            % to get list of item
submanipFieldName={'ratingReward', 'ratingRewardPicture', 'ratingPunishment', 'ratingEffort'}; % to use as field name in data
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
        subdata.condition.trialNumber     = allData(:, 1)';
        subdata.condition.itemNumber      = allData(:, 2)';                                      % item number
        subdata.condition.itemSubtype     = allData(:, 2)'>12;                                   % item type
        
        
        
        %% BEHAVIOUR
        % ===========================================================================
        
        subdata.behavior.rating  = allData(:, 3)';                  % maximal force recorded during the trial (in N or abstract unit, depending on the device)
        
        %% MISC
        % ===========================================================================
        for iSession = 1:nSession
            subdata.misc.session(iSession).fileName = fileList(iSession).fileName ;
        end
        if isfield(fileList(iSession), typeItemList{iSubmanip})
            subdata.misc.itemList =  getfield(fileList(iSession), typeItemList{iSubmanip});
        end
        
    catch err
        subdata = [];
    end
    data.(submanipFieldName{iSubmanip}) = subdata;

end


end