function [dynamo,behavior] = process_dynamo_accu(gripdata,subdata)
% process_dynamo - 
%
%
% Inputs:
%    gripdata - 
%
% Outputs:
%   dynamo - 
%   behavior - 
%

% parameters
%%% filtering
a = [0.5 0.5]; % transfer function
b = 1; 
%%% criterions
forcePeakCriterion = 0.10;      % exclude non participated trials ( forcePeak < 5% of maximal observed force)
reactionTimeCriterion = 0.20;   % detect reaction time ( first time marker with yank > 20% of yankPeak)
reactionTimeLowerBound = 0.200; % exclude reaction time estimates artefact( rt < 200 ms. )
onsetCriterion = 0.50;
%%% parameters
calib = subdata.calibration.forcePeak;


% add events 
ntrial = numel(gripdata.grip);
for i = 1:ntrial
    gripdata.nTrial{i} = i*ones(1,numel(gripdata.grip{i}));
    gripdata.responseOnset{i} = zeros(1,numel(gripdata.grip{i}));
    gripdata.responseOnset{i}(1) = 1;
end

% format into structure
dynamo.raw_force = [gripdata.grip{:}];
dynamo.time = [gripdata.time{:}] - gripdata.time{1}(1);
dynamo.nTrial = [gripdata.nTrial{:}];
dynamo.responseOnset = [gripdata.responseOnset{:}];
dynamo.force = nan(1,numel(dynamo.raw_force));
fmax = max(dynamo.raw_force);

% format into a table
behaviorNames = {'nTrial','forcePeak','norm_forcePeak','norm_forceSum','isCorrect',...
                  'time2peak','velocityPeak','norm_velocityPeak','relaxationPeak','norm_relaxationPeak',...
                  'responseTime',...
                  'forceDuration','preparatoryDuration',...
                  'trialDuration','norm_forceAmplitude','correct_force'};
nBehavior = numel(behaviorNames);
behavior = array2table(nan(ntrial,nBehavior),'VariableNames',behaviorNames);

targetForce = calib(subdata.condition.sessionNumber,:);
targetForce = targetForce(subdata.condition.handSide)';
targetForce = (subdata.condition.costValue').*targetForce;



for i=1:ntrial
    % extract
    raw = dynamo.raw_force(dynamo.nTrial==i);
    time = dynamo.time(dynamo.nTrial==i);
    
    % filter
    T = time(end) - time(1);
    dt = T/ numel(raw);
    f_force = filtfilt(a,b,raw);
    dynamo.force(dynamo.nTrial==i) = f_force;
    [fpeak,ipeak] = max(f_force);
    
    % extract
    behavior.nTrial(i) = i;
    behavior.forcePeak(i) = fpeak;
    behavior.norm_forcePeak(i) = fpeak/targetForce(i);
    behavior.norm_forceSum(i) = sum(f_force)/targetForce(i);
    behavior.isCorrect(i) = double(fpeak >= forcePeakCriterion*targetForce(i));

    if behavior.isCorrect(i)==1
        
        % dynamic metrics
        velocity = [0,diff(f_force)];
        velocity = filtfilt(a,b,velocity)/dt;
        dynamo.velocity(dynamo.nTrial==i) = velocity;
        [vpeak,jpeak] = max(velocity);
        rpeak = min(velocity);
        behavior.velocityPeak(i) = vpeak;
        behavior.norm_velocityPeak(i) = vpeak/fmax;
        behavior.relaxationPeak(i) = rpeak;
        behavior.norm_relaxationPeak(i) = rpeak/fmax;
        
        % estimate behavior times
        try
            try
                % tangent method
                acc = [0,diff(velocity)];
                acc = filtfilt(a,b,acc)/dt;

                % onset transition
                inflexion = (acc <=0 & f_force>= forcePeakCriterion*fpeak );
                j = find(inflexion,1,'first');
                beta = velocity(j);
                tangent = f_force(j) + beta*(time-time(j));
                baseline = f_force(1);
                irt = find(tangent>=baseline,1,'first');
                rt = time(irt)-time(1);

                % plateau transition
                plateau = fpeak;
                ionset = find(tangent>=plateau,1,'first');

                inflexion = (acc <=0 & velocity<=0 );
                j = find(inflexion,1,'last');
                beta = velocity(j);
                tangent = f_force(j) + beta*(time-time(j));
                ioffset = find(plateau>=tangent,1,'first');
            catch
                rt = nan;
            end
            if rt<=reactionTimeLowerBound
                rt = nan;
            end
            behavior.responseTime(i) = rt;
            t2p = time(ipeak) - rt - time(1);
            behavior.time2peak(i) = t2p;

            % estimate performance metrics
            onset = (f_force>= onsetCriterion*targetForce(i));
            ieff = find(onset,1,'first');
            behavior.forceDuration(i) = time(end) - time(ieff);
            behavior.preparatoryDuration(i) = time(ieff) - time(1);
            behavior.trialDuration(i) = time(end) - time(1);
            behavior.norm_forceAmplitude(i) = mean(f_force(ionset:ioffset)/targetForce(i));
            behavior.correct_force(i) = mean(f_force(ionset:ioffset)>=targetForce(i));
        end
    end
    
end