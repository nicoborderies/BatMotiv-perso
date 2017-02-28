
function [data] = load_data_mentalRP(subjectDir, sessionList)
% load_data_mentalRP extracts data for the manip mentalRP / mentalgripRP of the MBB Battery
% [data] = load_data_mentalRP(subjectDir, sessionList)
%
% IN
%       - subjectDir is the source directory
%       - sessionList is an optional argument to select session (by default, all sessions are extracted)
%           in "dir" order. ie. sessionList=[5 2] will extract data in the 5th and the 2d files listed in natural order (as session 1 and session2)
% OUT
%       - data containing the following field
%       -------------------------------------------------------------------------------------------------------------
%           * timeAllowedToAnswer      : Actual time allowed to the subject to answer
%       -------------------------------------------------------------------------------------------------------------
%           * condition: experimental condition
%               - sessionNumber        : converted to 1 to nSession
%               - trialNumber          : in a given session
%               - incentiveSign        : [1] trying to gain money / [-1] = trying to avoid losses
%               - incentiveLevel       : [1-6] Incentive level (i rank)
%               - incentiveValue       : Incentive value (in €)
%       -------------------------------------------------------------------------------------------------------------
%           * behavior: depending on subject choie
%               - nGoodAnswer          : number of good answer during the trial
%               - nTotalAnswer         : total number of anwers
%               - trialGain            : money earned during the trial
%               - totalGain            : total money earned since the beginning of the experiment
%               - rawAnswer            : Structure containing all answers recorded for each trial ([1] good answer / [0] bad answer)
%               - rawTimeValue         : Structure containing all answers onset times
%
%               Adjusted stats
%                   1/ because of time penalty, some answers were collected after calibTime : correct for this : nCorrectGoodAnswer & nCorrectedTotalAnswer
%               - nCorrectedGoodAnswer   : nGoodAnswer corrected for answers given after  timeAllowedToAnswer
%               - nCorrectedTotalAnswer: nTotalAnswer corrected for answers given after  timeAllowedToAnswer
%                   2/ Try to find a clever way to measure performance beyond nCorrectGoodAnswer
%               -nAdjustedGoodAnswer   : nCorrectedGoodAnswer + non used time (i.e. time after the last answer) divided by the mean time used per good answer (maximum adjustment = 1 )
%               -timePerGoodAnswer     : last answer onset time / nCorrectedGoodAnswer
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
MANIP_NAME = 'mental*sess'; % please keep the underscores for compatibility between versions !

% load files in a struct
if nargin<2, sessionList = []; end
[fileList, nSession] = tools.load_file(subjectDir, MANIP_NAME, sessionList);


%% Compute adjusted nGoodAnswer
% 1/ because of time penalty, some answers were collected after calibTime : correct for this : nCorrectGoodAnswer & nCorrectedTotalAnswer
% 2/ Try to find a clever way to measure performance beyond nCorrectGoodAnswer
% a/ nAdjustedGoodAnswer : add to nCorrectedGoodAnswer non used time (i.e. time after the last answer) divided by the mean time used per good answer (maximum adjustment = 1 )
% b/ timePerGoodAnswer : last answer time / nCorrectedGoodAnswer

for iSession = 1:nSession
    nTrial=max(fileList(iSession).data(:, 1));
    calibTime=fileList(iSession).calibtime;
    nGoodAnswer  = fileList(iSession).data(:, 3)';
    nTotalAnswer = fileList(iSession).data(:, 4)';
    rawAnswer=fileList(iSession).datatask.answer;
    rawAnswerTime=fileList(iSession).datatask.answertime;
    for iTrial=1:nTrial
        if nGoodAnswer(iTrial) == 0
            nCorrectedGoodAnswer(iTrial)=0;
            nAdjustedGoodAnswer(iTrial)=0;           
            nCorrectedTotalAnswer(iTrial)=nTotalAnswer(iTrial);
            timePerGoodAnswer(iTrial)=NaN;
        else
            trialRawAnswer=rawAnswer{iTrial};
            trialRawAnswerTime=rawAnswerTime{iTrial};
            lastAnswerTime=max(trialRawAnswerTime)-fileList(iSession).data(iTrial,9);
            unusedTime=calibTime-lastAnswerTime;
            if unusedTime < 0 % Answer after calibTime --> Remove the last answer
                nCorrectedTotalAnswer(iTrial)= nTotalAnswer(iTrial)-1; 
                if trialRawAnswer(end)==1
                    nCorrectedGoodAnswer(iTrial)=nGoodAnswer(iTrial)-1;
                else
                    nCorrectedGoodAnswer(iTrial)=nGoodAnswer(iTrial);
                end
                lastAnswerTime=trialRawAnswerTime(end-1)-fileList(iSession).data(iTrial,9);
                if trialRawAnswer(end-1)==0
                    unusedTime=calibTime-(lastAnswerTime+0.1*calibTime); % Add penalty time
                else
                    unusedTime=calibTime-lastAnswerTime;
                end
                unusedTime=max(unusedTime,0);
            else
                nCorrectedGoodAnswer(iTrial)=nGoodAnswer(iTrial);
                nCorrectedTotalAnswer(iTrial)=nTotalAnswer(iTrial);
            end
            if nCorrectedGoodAnswer(iTrial)==0
                timePerGoodAnswer(iTrial)= NaN;
                nAdjustedGoodAnswer(iTrial)=0; 
            else
                timePerGoodAnswer(iTrial)= lastAnswerTime / nCorrectedGoodAnswer(iTrial);
                nAdjustedGoodAnswer(iTrial) = nCorrectedGoodAnswer(iTrial) + min(unusedTime / timePerGoodAnswer(iTrial),1);
            end
        end
    end
    fileList(iSession).data=[fileList(iSession).data nCorrectedGoodAnswer' nCorrectedTotalAnswer' nAdjustedGoodAnswer' timePerGoodAnswer'];
end
    

% concatenate sessions
allData = vertcat(fileList.data);


%% Actual time allowed to the subject to answer
% ===========================================================================

for iSession = 1:nSession
    data.timeAllowedToAnswer(iSession)=fileList.calibtime;   % Actual time allowed to the subject to answer
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

% 1/ because of time penalty, some answers were collected after calibTime : correct for this : nCorrectGoodAnswer & nCorrectedTotalAnswer
% 2/ Try to find a clever way to measure performance beyond nCorrectGoodAnswer
% a/ nAdjustedGoodAnswer : add to nCorrectedGoodAnswer non used time (i.e. time after the last answer) divided by the mean time used per good answer (maximum adjustment = 1 )
% b/ timePerGoodAnswer : last answer time / nCorrectedGoodAnswer

data.behavior.nGoodAnswer  = allData(:, 3)';                % number of good answer during the trial
data.behavior.nTotalAnswer  = allData(:, 4)';               % total number of anwers
data.behavior.nCorrectedGoodAnswer  = allData(:, 12)';      % corrected number of good answer during the trial (remove answer after calib time)
data.behavior.nCorrectedTotalAnswer  = allData(:, 13)';     % corrected total number of anwers (remove answer after calib time)
data.behavior.nAdjustedGoodAnswer  = allData(:, 14)';       % Adjusted number of good answer during the trial
data.behavior.timePerGoodAnswer  = allData(:, 15)';         % Last answer time / nCorrectedGoodAnswer
data.behavior.trialGain  = allData(:, 6)';                  % money earned during the trial
data.behavior.totalGain  = allData(:, 7)';                  % total money earned since the beginning of the experiment

data.behavior.rawAnswer={};
data.behavior.rawAnswerTime={};

for iSession = 1:nSession
    data.behavior.rawAnswer={data.behavior.rawAnswer{1:end} fileList(iSession).datatask.answer{1:end}};
    data.behavior.rawAnswerTime={data.behavior.rawAnswerTime{1:end} fileList(iSession).datatask.answertime{1:end}};
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



