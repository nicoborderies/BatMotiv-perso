function [result] = process_data_discount_modelfree(dataStructure, option)
% process_data_discount_modelfree do the model free analysis for the task discount of the MBB Battery for a single subject
% [data] = process_data_discount_modelfree(dataStructure)
%
% IN
%       - data structure from the load_data_discout function
%       - option.nBin = number of bin for  singleStat 
%
% OUT
%       - results containing one field by submanip. Each submanip contains:
%       -------------------------------------------------------------------------------------------------------------
%           * singleStat : statistics usefull for a single subject (e.g. to plot) but won't be passed to the second level... ex variance, distribution etc
%               - binOfDifferenceValue               : mean value of each bin (of difference rating between the immediate and the delayed item) (1 * 6 bins * nSessions)
%               - binOfDifferenceValueDelayedChoice  : percentage of choice of the delayed item by bin (of difference rating between the immediate and the delayed item) (1 * 6 bins * nSession)
%               - binOfDelay                         : mean value of each bin (of delay)  (1* 6 bins*nSession )
%               - binOfDelayDelayedChoice            : percentage of choice of the delayed item by bin (of delay) (1* 6 bins*nSession )
%           * firstLevelStat : statistics easy to pass to the second level (e.g for t.tests)
%               - meanDelayedChoice                  : % of choice of the delayed item (1*1*nSession)
%               - betaDifferenceValue                : beta of the logistic regression (1*1*nSession) : difference value regressor
%               - betaDelay                          : beta of the logistic regression (1*1*nSession) : delay regressor
%               * Just for coherence with weighting functions
%               - betaDifferenceValueFromSeparateGLM : beta of the logistic regression (1*1*nSession) : difference value regressor
%               - betaDelayFromSeparateGLM           : beta of the logistic regression (1*1*nSession) : delay regressor

if nargin<2; option.nBin=6; end
submanipFieldName = fieldnames(dataStructure);
nSubManip = length(submanipFieldName);
result=struct;

selectSubManip=[];
for iSubManip = 1:nSubManip
    if ~isempty(dataStructure.(submanipFieldName{iSubManip}))
        selectSubManip=[selectSubManip iSubManip];
    end
end

%% ANALYSES
% ===========================================================================
% Loop across submanip
for iSubManip = selectSubManip
    % subset to submanip data
    subData=dataStructure.(submanipFieldName{iSubManip});
    % rename for clarity
    differenceValue = subData.condition.differenceRatingValue;
    constant=ones(size(differenceValue));
    delay = subData.condition.delay;
    isDelayedItemChoice = subData.behavior.isDelayedItemChoice;
    differenceValueBin=[];
    delayBin=[];
    for iSession = unique(subData.condition.sessionNumber)
        nTrial= sum(subData.condition.sessionNumber==iSession);
        % select data from the current session
        sessionDifferenceValue = differenceValue(subData.condition.sessionNumber==iSession);
        sessionDelay = delay(subData.condition.sessionNumber==iSession);
        sessionIsDelayedItemChoice = isDelayedItemChoice(subData.condition.sessionNumber==iSession);
        sessionDifferenceValueBin=ones(1,nTrial);
        sessionDelayBin=ones(1,nTrial);
        % make bin
        quantileDifferenceValueBin=quantile(sessionDifferenceValue, (1/option.nBin):(1/option.nBin):((option.nBin-1)/option.nBin));
        quantileDelay=quantile(sessionDelay, (1/option.nBin):(1/option.nBin):((option.nBin-1)/option.nBin));
        for iQuantile=1:length(quantileDifferenceValueBin)
             sessionDifferenceValueBin(sessionDifferenceValue>quantileDifferenceValueBin(iQuantile))=sessionDifferenceValueBin(sessionDifferenceValue>quantileDifferenceValueBin(iQuantile))+1;
             sessionDelayBin(sessionDelay>quantileDelay(iQuantile))=sessionDelayBin(sessionDelay>quantileDelay(iQuantile))+1;
        end
        differenceValueBin=[differenceValueBin sessionDifferenceValueBin];
        delayBin=[delayBin sessionDelayBin];
     
        % do linear fit
        regressor=glmfit([zscore(sessionDifferenceValue)' zscore(sessionDelay)'],sessionIsDelayedItemChoice', 'binomial', 'link', 'logit');
        subResult.firstLevelStat.betaDifferenceValue(1,1,iSession)=regressor(2);
        subResult.firstLevelStat.betaDelay(1,1,iSession)=regressor(3);
        regressor=glmfit([zscore(sessionDifferenceValue)'],sessionIsDelayedItemChoice', 'binomial', 'link', 'logit');
        subResult.firstLevelStat.betaDifferenceValueFromSeparateGLM(1,1,iSession)=regressor(2);
        regressor=glmfit([zscore(sessionDelay)'],sessionIsDelayedItemChoice', 'binomial', 'link', 'logit');
        subResult.firstLevelStat.betaDelayFromSeparateGLM(1,1,iSession)=regressor(2);
    end
    % fill subresult structure
    subResult.singleStat.binOfDifferenceValue  = tools.tapply(differenceValue, {constant, differenceValueBin,subData.condition.sessionNumber}, @mean);
    subResult.singleStat.binOfDifferenceValueDelayedChoice = tools.tapply(isDelayedItemChoice, {constant, differenceValueBin,subData.condition.sessionNumber}, @mean);
    subResult.singleStat.binOfDelay = tools.tapply(delay, {constant, delayBin, subData.condition.sessionNumber}, @mean);
    subResult.singleStat.binOfDelayDelayedChoice = tools.tapply(isDelayedItemChoice, {constant, delayBin,subData.condition.sessionNumber}, @mean);
    subResult.firstLevelStat.meanDelayedChoice = tools.tapply(isDelayedItemChoice, {constant,constant,subData.condition.sessionNumber}, @mean);
    % fill result structure
    result.(submanipFieldName{iSubManip})= subResult;
end



