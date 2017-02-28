
function [data] = load_data_gripR(subjectDir, sessionList)
% LOAD_DATA_gripR extracts data for the manip gripR of the MBB Battery
% [data] = load_data_gripR(subjectDir, sessionList)
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
%               - incentiveSign        : [1] trying to gain money / [-1] = trying to avoid losses
%               - incentiveLevel       : [1-6] Incentive level (i rank)
%               - incentiveValue       : Incentive value (in €)
%       -------------------------------------------------------------------------------------------------------------
%           * behavior: depending on subject choie
%               - forcePeak            : maximal force recorded during the trial (in N or abstract unit, depending on the device)
%               - forceSum             : sum of forces recorded during the trial (in N or abstract unit, depending on the device)
%               - normalizedForcePeak  : normalized maximal force  (in % of calibration force)
%               - trialGain            : money earned during the trial
%               - totalGain            : total money earned since the beginning of the experiment
%               - rawForceValue        : Structure containing all force values recorded for each trial
%               - rawTimeValue         : Structure containing all timing values recorded for each trial
%       -------------------------------------------------------------------------------------------------------------
%           * calibration:
%               - forcePeak            : Maximal force recorded during the calibration (1 value / session)
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
%                   - fileName         : name of the loaded file

%% LOAD FILES
% ===========================================================================
% result files should include manip 'manipName'
MANIP_NAME = '_gripR_s'; % please keep the underscores for compatibility between versions !

% load files in a struct
if nargin<2, sessionList = []; end
[fileList, nSession] = tools.load_file(subjectDir, MANIP_NAME, sessionList);

% concatenate sessions
allData = vertcat(fileList.data);


%% CALIBRATION
% ===========================================================================

for iSession = 1:nSession
    data.calibration.forcePeak(iSession)=fileList.calib;   % maximal force recorded during calibration
end

%% CONDITIONS
% ===========================================================================
data.condition.sessionNumber=[];
for iSession = 1:nSession
    [nTrial,~]=size(fileList(iSession).data);
    data.condition.sessionNumber =   [data.condition.sessionNumber iSession * ones(1,nTrial)] ;
end
data.condition.trialNumber     = allData(:, 1)';
data.condition.incentiveSign   = sign(allData(:, 2)');      % 1 = trying to gain money / -1 = trying to avoid losses
data.condition.incentiveLevel  = abs(allData(:, 2)');       % Incentive level (rank)
incentiveList=[0.01 0.2 0.5 1 5 20];
data.condition.incentiveValue  = incentiveList(data.condition.incentiveLevel);       % Incentive value (in €)

%% BEHAVIOUR
% ===========================================================================
                                        
data.behavior.forcePeak  = allData(:, 3)';                  % maximal force recorded during the trial (in N or abstract unit, depending on the device)
data.behavior.forceSum  = allData(:, 4)';                   % sum of forces recorded during the trial (in N or abstract unit, depending on the device)
data.behavior.normalizedForcePeak  = allData(:, 5)';        % normalized maximal force  (in % of calibration force)
data.behavior.trialGain  = allData(:, 6)';                  % money earned during the trial
data.behavior.totalGain  = allData(:, 7)';                  % total money earned since the beginning of the experiment

data.behavior.rawForceValue={};
data.behavior.rawTimeValue={};

for iSession = 1:nSession
    data.behavior.rawForceValue={data.behavior.rawForceValue{1:end} fileList(1).gripdata.grip{1:end}};
    data.behavior.rawTimeValue={data.behavior.rawTimeValue{1:end} fileList(1).gripdata.time{1:end}};
end


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
for iSession = 1:nSession
    data.misc.session(iSession).fileName = fileList(iSession).fileName ;
end



