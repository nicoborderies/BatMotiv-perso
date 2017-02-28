
function [data] = load_data_weight(subjectDir, sessionList)
% LOAD_DATA_discount extracts data for all weight manip of the MBB Battery
% [data] = load_data_weight(subjectDir, sessionList)
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
%               - benefitItemNumber                        : benefit item number : reward in RP & RE, (avoid) punishment in PE
%               - costItemNumber                           : cost item number : effort in RE & PE, punishment in RP
%               - benefitItemRatingValue                   : rating of the benefit item number : reward in RP & RE, (avoid) punishment in PE
%               - costItemRatingValue                      : rating of the cost item number : effort in RE & PE, punishment in RP
%       -------------------------------------------------------------------------------------------------------------
%           * behavior: depending on subject choie
%               - RT                   : reaction time
%               - sideChoice           : choice side [1] right [-1] left
%               - isGoChoice           : % choice type: [1] go : accept the trade-off [0] no go : refuse the trade-off
%       -------------------------------------------------------------------------------------------------------------
%           * misc
%               - session
%                   - fileName         : name of the loaded file
%                   - benefitList : actual items used
%                   - costList    : actual items used


data=struct;

subManipList={'WeightRE_s', 'WeightPE_s', 'WeightRP_s'};                                    % to use as MANIP_NAME
typeItemList={'rewardlist', 'rewardlist' 'punishlist', 'effortlist'};                       % to get list of item
submanipFieldName={'weightRE', 'weightPE', 'weightRP'};                                     % to use as field name in data
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
        subdata.condition.benefitItemNumber          = allData(:, 2)';                         % benefit item number : reward in RP & RE, (avoid) punishment in PE
        subdata.condition.costItemNumber             = allData(:, 4)';                         % cost item number : effort in RE & PE, punishment in RP
        subdata.condition.benefitItemRatingValue     = allData(:, 3)';                         % rating of the benefit item number : reward in RP & RE, (avoid) punishment in PE
        subdata.condition.costItemRatingValue        = allData(:, 5)';                         % rating of the cost item number : reward in RP & RE, (avoid) punishment in PE
        
         
        %% BEHAVIOUR
        % ===========================================================================
        
        subdata.behavior.RT  = allData(:, 6)';                         % RT
        subdata.behavior.sideChoice  = allData(:, 7)';                 % choice side: [1] right [-1] left
        subdata.behavior.isGoChoice  = allData(:, 8)';                 % choice type: [1] go : accept the trade-off [0] no go : refuse the trade-off
        
        %% MISC
        % ===========================================================================
        for iSession = 1:nSession
            subdata.misc.session(iSession).fileName = fileList(iSession).fileName ;
        end
        iCompteur=0;
        for iFieldType = 1:4
            if isfield(fileList(iSession), typeItemList{iFieldType})
                if iCompteur == 0
                    subdata.misc.benefitList =  getfield(fileList(iSession), typeItemList{iFieldType});
                    iCompteur=iCompteur+1;
                else
                    subdata.misc.costList =  getfield(fileList(iSession), typeItemList{iFieldType});
                end
            end
        end

    catch 
        subdata=[];
    end
    
    data.(submanipFieldName{iSubmanip})=subdata;

    
end
end