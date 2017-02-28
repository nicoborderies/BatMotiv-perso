function [result] = process_learning_modelfree(subjectData, option) 
 % LOAD_DATA_learning extracts data for the manip learning of the MBB Battery 
 % [data] = process_learning_modelfree(subjectDir, sessionList) 
 % 
 % IN 
 %       - data structure from the load_data_mentalIAPS function 
 %       - option is an optionnal astructure with the following fields 
 %              * nBin : number of bin (rounded to the closest possibility) 
 % 
 % OUT 
 %       - results containing the following field 
 %       ------------------------------------------------------------------------------------------------------------- 
 %           * singleStat : statistics usefull for a single subject (e.g. to plot) but won't be passed to the second level... ex variance, distribution etc 
 %               * positiveLearningCurve : learning curve : a 1 * nBin * nSession matrix (positive cue) 
 %               * negativeLearningCurve : learning curve : a 1 * nBin * nSession matrix (negative cue) 
 % 
 %           * firstLevelStat : statistics easy to pass to the second level (e.g for t.tests) 
 %               *meanCorrectByValence : 1*2*nSession (first column: negative cue; second column : positive cue) 
 %               * betaLearning                   : from the logistic regression isOptimalChoice ~ trialNumber (in a given valence) * valence 
%               * betaValence                    : from the logistic regression isOptimalChoice ~ trialNumber (in a given valence) * valence 
 %               * betaInteractionLearningValence : from the logistic regression isOptimalChoice ~ trialNumber (in a given valence) * valence 
 %               * betaLearningPositive           : from the logistic regression isOptimalChoice ~ trialNumber (Positive pair only) 
 %               * betaLearningNegative           : from the logistic regression isOptimalChoice ~ trialNumber (Negative pair only) 
 % 
 % 
 % NB : sessionNumber is converted to get the real number of SESSIONS (i.e. days) 
 %    : blockNumber is the number of block by session 
 
 
 
 
 option = tools.check_option(option, {... 
     'nBin',6 ... 
     }); 
      
 % extract learning 
 [condition,behavior]=tools.extractManip(subjectData,{'learning'}); 
 
 
 
 
 % Get the real "session" number 
 sessionNumber = floor((condition.sessionNumber-1)/max(condition.blockNumber))+1; 
 nSession=max(sessionNumber); 
 nBlock=max(condition.blockNumber); 
 
 
 %Rename and subset to non-neutral pair 
 sessionNumber   = sessionNumber(condition.pairValence~=0); 
 blockNumber     = condition.blockNumber(condition.pairValence~=0); 
 isOptimalChoice = behavior.isOptimalChoice(condition.pairValence~=0); 
 pairValence     = condition.pairValence(condition.pairValence~=0); 
 constant=ones(size(sessionNumber)); 
 
 
 %Compute bin size 
 bin=zeros(size(sessionNumber)); 
 nTrialByValence=length(sessionNumber)/(nSession*nBlock*2); 
 binSize = round(nTrialByValence/option.nBin); 
 nBin=nTrialByValence/binSize; 
 bin(pairValence==-1)=repmat(sort(repmat(1:nBin, 1, binSize)),1, nSession*nBlock); 
 bin(pairValence==1)=repmat(sort(repmat(1:nBin, 1, binSize)),1, nSession*nBlock); 
 
 
 %% ANALYSES 
 % =========================================================================== 
 
 
 %% do singleStat 
 learningCurve=tools.tapply(isOptimalChoice, {pairValence, bin, sessionNumber}, @mean); 
 result.singleStat.negativeLearningCurve = learningCurve(1,:,:); 
 result.singleStat.positiveLearningCurve = learningCurve(2,:,:); 
 %% do firstLevelStat 
 bin(pairValence==-1)=repmat(1:nTrialByValence,1, nSession*nBlock); 
 bin(pairValence==1) =repmat(1:nTrialByValence,1, nSession*nBlock); 
 result.firstLevelStat.meanCorrectByValence = tools.tapply(isOptimalChoice, {constant, pairValence,sessionNumber}, @mean); 
 for iSession = unique(sessionNumber) % Loop across session 
    % select data from the current session 
    sessionPairValence=pairValence(sessionNumber==iSession); 
    sessionBin=bin(sessionNumber==iSession); 
    sessionIsOptimalChoice=isOptimalChoice(sessionNumber==iSession); 
    %fit 
    regressor=glmfit([sessionBin; sessionPairValence; sessionBin.*sessionPairValence]',sessionIsOptimalChoice', 'binomial', 'link', 'logit'); 
     result.firstLevelStat.betaLearning(1,1,iSession) = regressor(2); 
     result.firstLevelStat.betaValence(1,1,iSession) = regressor(3); 
     result.firstLevelStat.betaInteractionLearningValence(1,1,iSession) = regressor(4); 
     %fit 
     regressor=glmfit([sessionBin(sessionPairValence==1)]',sessionIsOptimalChoice(sessionPairValence==1)', 'binomial', 'link', 'logit'); 
     result.firstLevelStat.betaLearningPositive(1,1,iSession) = regressor(2); 
     regressor=glmfit([sessionBin(sessionPairValence==-1)]',sessionIsOptimalChoice(sessionPairValence==-1)', 'binomial', 'link', 'logit'); 
     result.firstLevelStat.betaLearningNegative(1,1,iSession) = regressor(2); 
 end 
 
 

