function [result] = process_mentalRP_modelfree(subjectData, option)
% process_data_mentalRP_modelfree do the model free analysis for the task mentalRP of the MBB Battery for a single subject
% [data] = process_mentalRP_modelfree(dataStructure)
%
% IN
%       - data structure from the load_data_mentalRP function
%       - option is an optionnal structure with the following fields
%              * incentiveMode : define how is express incitation 
%                      {'rank'}    = incentive is defined by its rank [1 to n levels] (default)
%                      {'money'}   = incentive is defined by actual money the suebject can win
%                      {@myFunction} = incentive is defined by myFunction(actual money)
%
% OUT
%       - results containing the following field
%       -------------------------------------------------------------------------------------------------------------
%           * singleStat : statistics usefull for a single subject (e.g. to plot) but won't be passed to the second level... ex variance, distribution etc
%               *[fitted variable]_By[factor] : is a 1*nLevel*nSession matrix where each cell is the mean value of the fitted variable for given session and level of the factor
%                fitted variables are :
%                     -nCorrectedGoodAnswer    = n of good answers corrected for answers given after  time allowed to answer
%                     -nCorrectedTotalAnswer   = total n of answers corrected for answers given after  time allowed to answer
%                     -nAdjustedGoodAnswer     = nCorrectedGoodAnswer + non used time (i.e. time after the last answer) divided by the mean time used per good answer (maximum adjustment = 1 )
%                     -timePerGoodAnswer       = last answer onset time / nCorrectedGoodAnswer
%                factor are
%                     - IncentiveLevel         = level of incentive (according to chosen option) whatever the sign (force to gain or to avoid loss)
%                     - Valence                = type of incentive (positive or negative)
%                     - PositiveIncentiveLevel = by level of incentive (positive incitation only)
%                     - NegativeIncentiveLevel = by level of incentive (positive incitation only)
%
%           * firstLevelStat : statistics easy to pass to the second level (e.g for t.tests)
%               *beta[factor]From_[fitted variable] : is the beta from the linear regression : fitted variable ~ factor
%                     -nCorrectedGoodAnswer    = n of good answers corrected for answers given after  time allowed to answer
%                     -nCorrectedTotalAnswer   = total n of answers corrected for answers given after  time allowed to answer
%                     -nAdjustedGoodAnswer     = nCorrectedGoodAnswer + non used time (i.e. time after the last answer) divided by the mean time used per good answer (maximum adjustment = 1 )
%                     -timePerGoodAnswer       = last answer onset time / nCorrectedGoodAnswer
%                factor are     
%                      - incentiveLevel         = level of incentive (according to chosen option) from GLM : y ~ Incentive level + valence
%                      - valence                = valence (according to chosen option) from GLM : y ~ Incentive level + valence
%                      - positiveIncentiveLevel = level of incentive restricted to valence ==  1 (try to gain)
%                      - negativeIncentiveLevel = level of incentive restricted to valence == -1 (try to avoid losses)
%
%
% NB : for regressions, regressors and fitted data are z-scored

option = tools.check_option(option, {...
    'incentiveMode',{'rank'} ...
    });
    
% extract learning
[condition,behavior]=tools.extractManip(subjectData,{'mentalRP'});

%% ANALYSES
  
% apply the same analysis to all these variables
listYFieldNames = {'nCorrectedGoodAnswer', 'nCorrectedTotalAnswer', 'nAdjustedGoodAnswer', 'timePerGoodAnswer'};
nY= length(listYFieldNames);


switch class(option.incentiveMode{1})
    case 'function_handle'
        myFunction=option.incentiveMode{1};
        incentiveLevel=myFunction(condition.incentiveValue);
    case 'char'
        switch option.incentiveMode{1}
            case 'rank'
                 incentiveLevel=condition.incentiveLevel;
            case 'money'
                 incentiveLevel=condition.incentiveValue;                
            otherwise
                warning('unknown incentive known')
                incentiveLevel=condition.incentiveLevel;
        end
    otherwise
        warning('unknown incentive known')
        incentiveLevel=condition.incentiveLevel;
end

% rename for clarity
valence=condition.incentiveSign;
sessionNumber=condition.sessionNumber;
constant=ones(size(sessionNumber));
for iY = 1:nY % loop acros variables to fit
    y = behavior.(listYFieldNames{iY});
    %% do single Stat
    result.singleStat.([listYFieldNames{iY} '_ByIncentiveLevel']) = tools.tapply(y, {constant, incentiveLevel, sessionNumber}, @nanmean);
    result.singleStat.([listYFieldNames{iY} '_ByValence']) = tools.tapply(y, {constant, valence,sessionNumber}, @nanmean);
    result.singleStat.([listYFieldNames{iY} '_ByPositiveIncentiveLevel']) = tools.tapply(y(valence==1), {constant(valence==1),incentiveLevel(valence==1),sessionNumber(valence==1)}, @nanmean);
    result.singleStat.([listYFieldNames{iY} '_ByNegativeIncentiveLevel']) = tools.tapply(y(valence==-1), {constant(valence==-1),incentiveLevel(valence==-1),sessionNumber(valence==-1)}, @nanmean);    
    %% do first level Stat 
    for iSession = unique(sessionNumber) % loop across sessions
        % select relevant part of the data
        sessionY=y(sessionNumber==iSession);
        sessionIncentiveLevel=incentiveLevel(sessionNumber==iSession);
        sessionValence=valence(sessionNumber==iSession);
        betaGLM1(1,:,iSession)=glmfit([zscore(sessionIncentiveLevel)' zscore(sessionValence)'],tools.nanzscore(sessionY)');
        betaGLM2(1,:,iSession)=glmfit([zscore(sessionIncentiveLevel(sessionValence==1))'],tools.nanzscore(sessionY(sessionValence==1))');
        betaGLM3(1,:,iSession)=glmfit([zscore(sessionIncentiveLevel(sessionValence==-1))'],tools.nanzscore(sessionY(sessionValence==-1))');
    end
    result.firstLevelStat.(['betaIncentiveLevelFrom_' listYFieldNames{iY}])= betaGLM1(:,2,:);
    result.firstLevelStat.(['betaValenceLevelFrom_' listYFieldNames{iY}])= betaGLM1(:,3,:);
    result.firstLevelStat.(['betaPositiveIncentiveLevelFrom_' listYFieldNames{iY}])= betaGLM2(:,2,:);
    result.firstLevelStat.(['betaNegativeIncentiveLevelFrom_' listYFieldNames{iY}])= betaGLM3(:,2,:);
end

end




