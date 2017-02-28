function [result] = process_data_choice_modelfree(dataStructure, option)
% process_data_choice_modelfree do the model free analysis for the task choice of the MBB Battery for a single subject
% [data] = process_data_choice_modelfree(dataStructure)
%
% IN
%       - data structure from the load_data_choice function
%       - option.nBin = number of bin for  singleStat 
%
% OUT
%       - results containing one field by submanip. Each submanip contains:
%       -------------------------------------------------------------------------------------------------------------
%           * singleStat : statistics usefull for a single subject (e.g. to plot) but won't be passed to the second level... ex variance, distribution etc
%               - binOfDifferenceRatingValue  : mean diffence rating values (1 * 6 bins *nSessions)
%               - binSideChoiceByValue        : percentage of "right-hand" choice by difference rating value (nSession * 6 bins)
%           * firstLevelStat : statistics easy to pass to the second level (e.g for t.tests)
%               - correctlyPredicted          : percentage of choices accuratly predicted by rating (1*1*nSession)
%               - betaDifferenceValue         : beta of the logistic regression (1*1*nSession) - (difference value is zscored)


if nargin<2; option.nBin=2; end
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
    differenceValue=subData.condition.differenceRatingValue;
    constant=ones(size(differenceValue));
    sideChoice=subData.behavior.sideChoice>0;
    if isempty(strfind(submanipFieldName{iSubManip}, 'eward')) % Not mandatory, it's just to plot all submanip in the same direction (choice increase with rating)
          %sideChoice=1-sideChoice;                             
          differenceValue = - differenceValue;                 % Change the sign of rating could be more intuitive ?
    end
    differenceValueBin=[];
    for iSession = unique(subData.condition.sessionNumber)
        nTrial= sum(subData.condition.sessionNumber==iSession);
        % select data from the current session
        sessionDifferenceValue=differenceValue(subData.condition.sessionNumber==iSession);
        sessionSideChoice=sideChoice(subData.condition.sessionNumber==iSession)>0;
        sessionDifferenceValueBin=ones(1,nTrial);
        % make bin
        quantileDifferenceValue=quantile(sessionDifferenceValue, (1/option.nBin):(1/option.nBin):((option.nBin-1)/option.nBin));
        for iQuantile=1:length(quantileDifferenceValue)
             sessionDifferenceValueBin(sessionDifferenceValue>=quantileDifferenceValue(iQuantile))=sessionDifferenceValueBin(sessionDifferenceValue>=quantileDifferenceValue(iQuantile))+1;
        end
        differenceValueBin=[differenceValueBin sessionDifferenceValueBin];
        % do linear fit
        regressor=glmfit(zscore(sessionDifferenceValue)',sessionSideChoice', 'binomial', 'link', 'logit');
        subResult.firstLevelStat.betaDifferenceValue(1,1,iSession)=regressor(2);
    end
    % fill subresult structure
    subResult.singleStat.binOfDifferenceRatingValue = tools.tapply(differenceValue, {constant, differenceValueBin, subData.condition.sessionNumber}, @mean);
    subResult.singleStat.binSideChoiceByValue = tools.tapply(sideChoice, {constant, differenceValueBin,subData.condition.sessionNumber}, @mean);
    predictedSide = differenceValue(differenceValue~=0)>0;
    subResult.firstLevelStat.correctlyPredicted = tools.tapply(predictedSide == sideChoice(differenceValue~=0), {constant(differenceValue~=0),constant(differenceValue~=0), subData.condition.sessionNumber(differenceValue~=0)}, @mean);
    % fill result structure
    result.(submanipFieldName{iSubManip})= subResult;
end



