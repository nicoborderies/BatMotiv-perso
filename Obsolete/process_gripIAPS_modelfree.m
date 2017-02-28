function [result] = process_gripIAPS_modelfree(subjectData, option)
% LOAD_DATA_gripIAPS extracts data for the manip gripIAPS of the MBB Battery
% [data] = load_data_gripIAPS(subjectDir, sessionList)
%
% IN
%       - data structure from the load_data_mentalIAPS function
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
%                     - IncentiveLevel         = level of incentive (according to chosen option)
%                     - pictureValence         = valence of the IAPS picture
%
%           * firstLevelStat : statistics easy to pass to the second level (e.g for t.tests)
%               *beta[factor]From_[fitted variable] : is the beta from the linear regression : fitted variable ~ IncentiveLevel + pictureArousal
%                     -normalizedForcePeak  : normalized maximal force  (in % of calibration force)
%                     -normalizedForceSum   : sum of normalized forces  (in % of calibration force)
%                     - zNormalizedForcePeak : zscored normalizedForcePeak
%                     - zNormalizedForceSum  : zscored normalizedForceSum
%                factor are     
%                     - IncentiveLevel         = level of incentive (according to chosen option)
%                     - pictureArousal         = absolute value of IAPS picture ([1] for emotionnal, [0] for neutral)
%
%
% NB : for regressions, regressors and fitted data are z-scored

option = tools.check_option(option, {...
    'incentiveMode',{'rank'} ...
    });
    
% extract learning
[condition,behavior]=tools.extractManip(subjectData,{'gripIAPS'});


%% ANALYSES

% apply the same analysis to all these variables
listYFieldNames = {'normalizedForcePeak', 'normalizedForceSum','zNormalizedForcePeak', 'zNormalizedForceSum'};
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
sessionNumber=condition.sessionNumber;
pictureValence=condition.pictureValence;
pictureArousal=abs(condition.pictureValence);
constant=ones(size(sessionNumber));
for iY = 1:nY % loop acros variables to fit
    y = behavior.(listYFieldNames{iY});
    %% do single Stat
    result.singleStat.([listYFieldNames{iY} '_ByIncentiveLevel']) = tools.tapply(y, {constant, incentiveLevel,sessionNumber}, @mean);
    result.singleStat.([listYFieldNames{iY} '_ByPictureValence']) = tools.tapply(y, {constant, pictureValence,sessionNumber}, @mean);  
    %% do first level Stat 
    for iSession = unique(sessionNumber) % loop across sessions
        % select relevant part of the data
        sessionY=y(sessionNumber==iSession);
        sessionIncentiveLevel=incentiveLevel(sessionNumber==iSession);
        sessionPictureArousal=pictureArousal(sessionNumber==iSession);
        betaGLM1(1,:, iSession)=glmfit([zscore(sessionIncentiveLevel)' zscore(sessionPictureArousal)'],zscore(sessionY)');
    end
    result.firstLevelStat.(['betaIncentiveLevelFrom_' listYFieldNames{iY}])= betaGLM1(:,2,:);
    result.firstLevelStat.(['betapictureArousalFrom_' listYFieldNames{iY}])= betaGLM1(:,3,:);
end
    

end




