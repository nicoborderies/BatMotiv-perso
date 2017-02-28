
function [data] = load_data_gripIAPS(subjectDir, sessionList)
% LOAD_DATA_gripIAPS extracts data for the manip gripIAPS of the MBB Battery
% [data] = load_data_gripIAPS(subjectDir, sessionList)
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
%               - pictureValence       : valence of the IAPS picture [-1] negative / [0] neutral / [1] positive
%       -------------------------------------------------------------------------------------------------------------
%           * behavior: depending on subject choie
%               - forcePeak            : maximal force recorded during the trial (in N or abstract unit, depending on the device)
%               - forceSum             : sum of forces recorded during the trial (in N or abstract unit, depending on the device)
%               - normalizedForcePeak  : normalized maximal force  (in % of calibration force)
%               - normalizedForceSum   : sum of normalized forces  (in % of calibration force)
%               - zNormalizedForcePeak : zscored normalizedForcePeak
%               - zNormalizedForceSum  : zscored normalizedForceSum
%               - trialGain            : money earned during the trial
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

%% INITIALIZE 
    % variable definition
    data=struct;
    data.interDimensions = struct;

    % session selection
%     if nargin<3;
        sessionList = []; 
%     end

    %% FILE EXTRACTION
    % ===========================================================================
    MANIP_NAME = 'GripIAPS'; % result files should include manip 'manipName'
    
    % test for file accessibility
     fileList=dir([subjectDir filesep '*' MANIP_NAME '*.mat']);
    if isempty(fileList)
        warning(['No file found for ' MANIP_NAME 'in directory: ' subjectDir])
        isFileFound =0;
    else
        isFileFound =1;
    end
    
        
    if isFileFound == 1 % check file existence
        [fileList, nSession] = tools.load_file(subjectDir, MANIP_NAME);

        
        %% MISC
        % ===========================================================================
        for iSession = 1:nSession
            data.misc.session(iSession).fileName = fileList(iSession).fileName ;
        end
        

        % compute morphological predicted Fmax
%         if ~isnan(option.morphology.antSkinFold)
%             subdata.misc.maxMorphologicalForce = Emax_morpho(option.morphology.antSkinFold , option.morphology.antSkinFold , option.morphology.foreArmCircumference , option.morphology.foreArmLength );
%         else
            data.misc.maxMorphologicalForce = NaN;
%         end
        
        % force signal baseline correction
        for iSession = 1:nSession
            [  data.misc.baseline , fileList(iSession).gripdata.grip ] = grip_correctbaseline( fileList(iSession).gripdata );
            for iT = 1:numel(fileList(iSession).data(:,3))
                fileList(iSession).data(iT,5) = max(fileList(iSession).gripdata.grip{1,iT});
                fileList(iSession).data(iT,6) = sum(safepos(fileList(iSession).gripdata.grip{1,iT})); % integrate force signal only for positive values
            end
        end
        
        % concatenate sessions
        allData = vertcat(fileList.data);
        
        
        %% CALIBRATION
        % ===========================================================================

        for iSession = 1:nSession
            data.calibration.forcePeak(iSession)=fileList(iSession).calib;   % maximal force recorded during calibration
        end

        data.calibration.maxObservedForcePeak= max([ data.calibration.forcePeak(iSession) , allData(:, 5)' ]);   % maximal force recorded during the whole task


        %% CONDITIONS
        % ===========================================================================
        subdata.condition.sessionNumber=[];
        for iSession = 1:nSession
            [nTrial,~]=size(fileList(iSession).data);
            subdata.condition.sessionNumber =   [subdata.condition.sessionNumber iSession * ones(1,nTrial)] ;
        end
            %  treatment
         if isfield(option.design,'subjectTable')
             ind = strfind(subjectDir,'sub');
             iSub =  str2num(subjectDir(ind+3:end));
             if option.design.subjectTable.sessionPlacebo(iSub)==1
                      subdata.condition.treatment = subdata.condition.sessionNumber ;
             else
                      subdata.condition.treatment = [subdata.condition.sessionNumber(subdata.condition.sessionNumber==2) , subdata.condition.sessionNumber(subdata.condition.sessionNumber==1)] ;
             end
         else
                 subdata.condition.treatment = subdata.condition.sessionNumber ;             
         end
        subdata.condition.trialNumber     = allData(:, 1)';
        subdata.condition.incentiveSign   = sign(allData(:, 2)');      % 1 = trying to gain money / -1 = trying to avoid losses
        subdata.condition.incentiveLevel  = abs(allData(:, 2)');       % Incentive level (rank)
        incentiveList=[0.01 0.2 0.5 1 5 20];
        subdata.condition.incentiveValue  = incentiveList(subdata.condition.incentiveLevel);       % Incentive value (in €)
        data.condition.pictureValence = allData(:, 3)';                                      % valence of the IAPS picture


        %% BEHAVIOUR
        % ===========================================================================
        % initialize
        subdata.behavior.rawForceValue={};
        subdata.behavior.rawTimeValue={};
        subdata.behavior.smoothForceValue={};
        subdata.behavior.normalizedForceValue={};
        subdata.behavior.smoothYankValue={};
        subdata.behavior.normalizedYankValue={};
        subdata.behavior.alignOnRT_normalizedForceValue = nan(numel(allData(:, 3)),60);
        subdata.behavior.correctResponseTiming = ones(1,numel(allData(:, 3)));
        subdata.behavior.time2forcePeak      = nan(1,numel(allData(:, 3)));
        subdata.behavior.time2yankPeak       = nan(1,numel(allData(:, 3)));
        subdata.behavior.rt                  = nan(1,numel(allData(:, 3)));
        subdata.behavior.yankPeak            = nan(1,numel(allData(:, 3)));
        subdata.behavior.normalizedYankPeak  = nan(1,numel(allData(:, 3)));
        
        
        
        % quality control 1
            % data selection
            forcePeakCriterion = 0.10;      % exclude non participated trials ( forcePeak < 10% of maximal observed force)
            forcePeakExclusion = ( allData(:,3)<forcePeakCriterion*max(allData(:,5)) ) ; 
            allData( forcePeakExclusion , 5) = NaN; 
            subdata.behavior.correctResponseTiming(forcePeakExclusion) = 0;
            responseTimingCriterion = 0.25; % detect premature response trials ( initial force > 25% of  forcePeak)
            reactionTimeCriterion = 0.20;   % detect reaction time ( first time marker with yank > 20% of yankPeak)
            reactionTimeLowerBound = 0.200; % exclude reaction time estimates artefact( rt < 200 ms. )

            
        % data extraction
        subdata.behavior.forcePeak  = allData(:, 5)';                  % maximal force recorded during the trial (in N or abstract unit, depending on the device)
        subdata.behavior.forceSum  = allData(:, 6)';                   % sum of forces recorded during the trial (in N or abstract unit, depending on the device)
        % subdata.behavior.normalizedForcePeak  = allData(:, 5)';        % normalized maximal force  (in % of calibration force)
        subdata.behavior.normalizedForcePeak  = (allData(:, 5)')./subdata.calibration.maxObservedForcePeak;        % normalized maximal force  (in % of maximal observed force)
        subdata.behavior.trialGain  = allData(:, 8)';                  % money earned during the trial
%         subdata.behavior.totalGain  = allData(:, 7)';                  % total money earned since the beginning of the experiment
        for iSession = 1:nSession
            subdata.behavior.rawForceValue = {subdata.behavior.rawForceValue{1:end} , fileList(iSession).gripdata.grip{1:end}};
            subdata.behavior.rawTimeValue = {subdata.behavior.rawTimeValue{1:end}  , fileList(iSession).gripdata.time{1:end}};
        end
        subdata.behavior.normalizedForceSum = subdata.behavior.forceSum/subdata.calibration.maxObservedForcePeak*100;

        
        
        % time-series preprocessing: temporal smoothing / detect premature responses / extract temporal derivative
        window=2; 
        h=ones(window,1)/window;
        for i=1:length(subdata.behavior.rawForceValue)
                trialDuration = subdata.behavior.rawTimeValue{i}(end)-subdata.behavior.rawTimeValue{i}(1);
                dt=(trialDuration)/length(subdata.behavior.rawForceValue{i});
                Fdynamic_temp=filter(h,1,subdata.behavior.rawForceValue{i});
                Fdynamic=fliplr(filter(h,1,fliplr(Fdynamic_temp)));
                subdata.behavior.smoothForceValue{i} = Fdynamic;
                subdata.behavior.normalizedForceValue{i} = Fdynamic./subdata.calibration.maxObservedForcePeak;
                
            if subdata.behavior.correctResponseTiming(i)==1
                [subdata.behavior.forcePeak(i),indexFmax(i)] = max(subdata.behavior.smoothForceValue{i});
                subdata.behavior.normalizedForcePeak(i)  = (subdata.behavior.forcePeak(i))./subdata.calibration.maxObservedForcePeak;     
                subdata.behavior.time2forcePeak(i) = subdata.behavior.rawTimeValue{i}(indexFmax(i))-subdata.behavior.rawTimeValue{i}(1);
                
                speed=diff(Fdynamic);
                Fspeed_temp=filter(h,1,speed);
                Fspeed=fliplr(filter(h,1,fliplr(Fspeed_temp)))/dt;
                subdata.behavior.smoothYankValue{i} = Fspeed;
                subdata.behavior.normalizedYankValue{i} = Fspeed./subdata.calibration.maxObservedForcePeak;

                if subdata.behavior.rawForceValue{1,i}(1) > responseTimingCriterion*max(allData(i,5)) % exclude velocity computation for premature trials (criterion: forcePeak after 25ms. > 20% of force peak)
                    subdata.behavior.yankPeak(i)=NaN;
                    subdata.behavior.normalizedYankPeak(i)=NaN;
                    subdata.behavior.correctResponseTiming(i)=0;
                else
                    subdata.behavior.yankPeak(i)=max(Fspeed);
                    subdata.behavior.normalizedYankPeak(i)=max(Fspeed)./subdata.calibration.maxObservedForcePeak;
                    [~,indexYmax(i)] = max(Fspeed);
                    subdata.behavior.time2yankPeak(i) = subdata.behavior.rawTimeValue{i}(indexYmax(i))-subdata.behavior.rawTimeValue{i}(1);
                    
                    % test for correct rt estimate
                    iter = 0; subdata.behavior.rt(i)=0;
                    while subdata.behavior.rt(i) < reactionTimeLowerBound
                        iter = iter + 1;
                        index = find(Fspeed >= reactionTimeCriterion*max(Fspeed),iter,'first');
                        indexRT(i) = index(iter);
                        subdata.behavior.rt(i) = subdata.behavior.rawTimeValue{i}(indexRT(i))-subdata.behavior.rawTimeValue{i}(1);
                    end
                    subdata.behavior.alignOnRT_normalizedForceValue(i,1:numel(subdata.behavior.normalizedForceValue{i}(indexRT(i):indexFmax(i)))) ...
                                                                                 = subdata.behavior.normalizedForceValue{i}(indexRT(i):indexFmax(i));
                    
                end
            end
        end
        
        
        % quality control 2: visualisation of  data processing
%         f = figure; hold on; set(f,'Color',[1 1 1]);
%         for iT = 1:numel(allData(:, 3))
%                 clf('reset');
%                 if subdata.behavior.correctResponseTiming(iT)==1
%                     color = 'k';   title('correct timing','Color',color);
%                 else
%                     color = 'r';   title('incorrect timing','Color',color);
%                 end
%                 
%                 subplot(2,1,1);    hold on;            
%                 h(1) =  plot(subdata.behavior.rawForceValue{iT},color);
%                 h(2) =  plot(subdata.behavior.smoothForceValue{iT},[color '--']);
%                 yy= ylim;
%                 if subdata.behavior.correctResponseTiming(iT)==1
%                     h(3) = scatter(indexFmax(iT),subdata.behavior.forcePeak(iT),'filled','r');
%                     if ~isnan(subdata.behavior.yankPeak(iT))
%                         h(4) = plot([ indexRT(iT)  indexRT(iT)],[yy(1) yy(2)],'r');
%                     end
%                 end
%                 ylabel('force)'); legend('raw','smoothed','fmax','rt');
% 
%                 subplot(2,1,2);  hold on;  
%                 h(1) =  plot(subdata.behavior.smoothYankValue{iT},[color '--']);
%                 if ~isnan(subdata.behavior.yankPeak(iT))
%                     h(2) = scatter(indexYmax(iT),subdata.behavior.yankPeak(iT),'filled','r');
%                 end
%                 ylabel('yank)'); legend('smoothed','ymax');
%                     
%                 xx=xlim;yy=ylim;text(xx(2),yy(2),num2str(iT));
%                 xlabel('time');
%                 title('correct timing');
%                 hold off;
%                 pause
%         end
        
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

        
    else
        data=[];
        
    end % file existence testing

    
%% TABLE



end

