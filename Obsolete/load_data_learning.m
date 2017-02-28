function [data] = load_data_learning(subjectDir, sessionList)
% LOAD_DATA_LEARNING extracts data for the learning_ox task of the MBB Battery
% [data] = load_data_learning(subjectDir, sessionList)
%
% IN
%       - subjectDir is the source directory
%       - sessionList is an optional argument to select session (by default, all sessions are extracted)
%           in "dir" order. ie. sessionList=[5 2] will extract data in the 5th and the 2d files listed in natural order (as session 1 and session2)
% OUT
%       - data containing the following field
%       -------------------------------------------------------------------------------------------------------------
%           * condition: experimental condition
%               - sessionNumber        : converted to 1 to nSession
%               - trialNumber          : in a given session
%               - pairValence          : [-1] negative / [0] neutral / [1] positive
%               - goodSide             : [-1] good cue on the left / [1] good cue on the right
%               - isOutcomePredictable : [0] unlikely outcome (25%) / [1] likely outcome (75%) / [NaN] neutral pair
%       -------------------------------------------------------------------------------------------------------------
%           * behavior: depending on subject choie
%               - choice               : [-1] left press / [1] right press
%               - isOptimalChoice      : [0] non-optimal choice / [1] optimal choice / [NaN] for neutral pair
%               - RT                   : reaction time
%               - outcome              : actual outcome [-1] -10euros / [0] 0euro / [1] 10euros
%       -------------------------------------------------------------------------------------------------------------
%           * oxymeter
%               - session
%                   - time             : time elapsed since the beginning of the session
%                   - tag              : [1] beginning of the trial / [2] ouctcome
%                   - heartRate        : Heart Rate
%                   - spO2             : Saturation
%       -------------------------------------------------------------------------------------------------------------
%           * misc
%               - session
%                   - type             : pair type
%                   - goodCueFile      : file used for the good cue (ie. leading to 10euros for the positive pair / 0euros for the negative pair)
%                   - badCueFile       : file used for the bad cue  (ie. leading to 0euros for the positive pair / -10euros for the negative pair)

%% LOAD FILES
% ===========================================================================
% result files should include manip 'manipName'
MANIP_NAME = 'learning_s';

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
data.condition.blockNumber = allData(:, 1)'  ;
data.condition.trialNumber = allData(:, 2)'  ;
data.condition.pairValence = allData(:, 3)'-2;                              % 1 negative / 2 neutral / 3 positive --> -1 / 0 / 1
data.condition.goodSide    = allData(:, 4)'  ;                              % -1 = good cue on the left, 1 = good cue on the right
data.condition.isOutcomePredictable = +(allData(:, 5)' > 0);                % -1 = unlikely outcome (25%), % 1 = likely outcome (75%)
data.condition.isOutcomePredictable(data.condition.pairValence==0)=NaN;     % There is no optimal choice for neutral pair


%% BEHAVIOUR
% ===========================================================================
data.behavior.choice  = allData(:, 7)';                                     % -1 = left press, right press
data.behavior.isOptimalChoice = +(allData(:, 8)>0)';                        % also encoded in allData(:, 8)
data.behavior.isOptimalChoice(data.condition.pairValence==0)=NaN;           % There is no optimal choice for neutral cue
data.behavior.RT      = allData(:, 6)';                                     % in s
data.behavior.outcome = allData(:, 10)';

%% OXYMETER
% ===========================================================================
% loop across sessions
for iSession = 1:nSession
    if isfield(fileList(iSession), 'oxdata')
        if ~isempty(fileList(iSession).oxdata);
    data.oxymeter.session(iSession).time      =  fileList(iSession).oxdata(:,12)';
    data.oxymeter.session(iSession).tag       =  fileList(iSession).oxdata(:,13)';
    data.oxymeter.session(iSession).heartRate =  fileList(iSession).oxdata(:,10)';
    data.oxymeter.session(iSession).spO2      =  fileList(iSession).oxdata(:,11)';
        end
    end
end

%% MISC
% ===========================================================================
% loop across sessions
for iSession = 1:nSession
    % filenames
    data.misc.session(iSession).fileName = fileList(iSession).fileName ;
    % get files used for the cues
    pairType={'negative', 'neutral', 'positive'};
    if isfield(fileList(iSession), 'cue_perm')
    for iCueType=1:3
        if fileList(iSession).cue_perm(iCueType)==1
            lS={'A', 'B'};
        else
            lS={'B', 'A'};
        end
        % store
        data.misc.session(iSession).pair(iCueType).type=pairType{iCueType};
        data.misc.session(iSession).pair(iCueType).goodCueFile=['Stim' num2str(fileList(iSession).sess_perm(iSession)) num2str(fileList(iSession).val_perm(iCueType)) lS{1} '.bmp'];
        data.misc.session(iSession).pair(iCueType).badCueFile =['Stim' num2str(fileList(iSession).sess_perm(iSession)) num2str(fileList(iSession).val_perm(iCueType)) lS{2} '.bmp'];
    end
    end
        
end




