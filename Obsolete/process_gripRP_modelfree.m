function [result] = process_gripRP_modelfree(subjectData, option)
% LOAD_DATA_gripRP extracts data for the manip gripRP of the MBB Battery
% [data] = process_gripRP_modelfree(subjectDir, sessionList)
%
% IN
%       - data structure from the load_data_mentalRP function
%       - option is an optionnal astructure with the following fields
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
%                     -normalizedForcePeak  : normalized maximal force  (in % of calibration force)
%                     -normalizedForceSum   : sum of normalized forces  (in % of calibration force)
%                factor are
%                     - IncentiveLevel         = level of incentive (according to chosen option) whatever the sign (force to gain or to avoid loss)
%                     - Valence                = type of incentive (positive or negative)
%                     - PositiveIncentiveLevel = by level of incentive (positive incitation only)
%                     - NegativeIncentiveLevel = by level of incentive (positive incitation only)
%
%           * firstLevelStat : statistics easy to pass to the second level (e.g for t.tests)
%               *beta[factor]From_[fitted variable] : is the beta from the linear regression : fitted variable ~ factor
%                     -normalizedForcePeak  : normalized maximal force  (in % of calibration force)
%                     -normalizedForceSum   : sum of normalized forces  (in % of calibration force)
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
[condition,behavior]=tools.extractManip(subjectData,{'gripRP'});

%% ANALYSES

% apply the same analysis to all these variables
listYFieldNames = {'normalizedForcePeak', 'normalizedForceSum'};
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




