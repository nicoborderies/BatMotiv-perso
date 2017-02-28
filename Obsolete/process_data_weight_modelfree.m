function [result] = process_data_weight_modelfree(dataStructure, option)
% process_data_weight_modelfree do the model free analysis for the task weight of the MBB Battery for a single subject
% [data] = process_data_weight_modelfree(dataStructure)
%
% IN
%       - data structure from the load_data_weight function
%       - option.nBin = number of bin for  singleStat 
%
% OUT
%       - results containing one field by submanip. Each submanip contains:
%       -------------------------------------------------------------------------------------------------------------
%           * singleStat : statistics usefull for a single subject (e.g. to plot) but won't be passed to the second level... ex variance, distribution etc
%               - binOfBenefitRatingValue     : mean value of each bin (of benefit rating) (1*6 bins*nSession)
%               - binOfBenefitChoice          : percentage of Go choice by bin (of benefit rating) (1*6 bins*nSession)
%               - binOfCostRating             : mean value of each bin (of cost rating) (1*6 bins*nSession)
%               - binOfCostChoice             : percentage of Go choice by bin (of cost rating) (1*6 bins*nSession)
%           * firstLevelStat : statistics easy to pass to the second level (e.g for t.tests)
%               - meanGo                      : % of go (1*1*nSession)
%               - betaBenefit                 : beta of the logistic regression (1*1*nSession) : benefit regresor (reward in RP & RE, (avoid) punishment in PE)
%               - betaCost                    : beta of the logistic regression (1*1*nSession) : cost regressor (effort in RE & PE, punishment in RP)
%               * Benefit and cost seem strongly anticorrelated (~ -0.8)... beta from separate GLMs are more coherent with data
%               - betaBenefitFromSeparateGLM  : beta of the logistic regression (1*1*nSession) : benefit regresor(reward in RP & RE, (avoid) punishment in PE)
%               - betaCostFromSeparateGLM     : beta of the logistic regression (1*1*nSession) : cost regressor (reward in RP & RE, (avoid) punishment in PE)

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
    benefitRating = subData.condition.benefitItemRatingValue;
    costRating = subData.condition.costItemRatingValue;
    isGoChoice = subData.behavior.isGoChoice;
    constant=ones(size(isGoChoice));
    benefitBin=[];
    costBin=[];
    for iSession = unique(subData.condition.sessionNumber) % Loop across session
        nTrial= sum(subData.condition.sessionNumber==iSession);
        % select data from the current session
        sessionBenefitRating = benefitRating(subData.condition.sessionNumber==iSession);
        sessionCostRating = costRating(subData.condition.sessionNumber==iSession);
        sessionIsGoChoice = isGoChoice(subData.condition.sessionNumber==iSession);
        sessionBenefitBin=ones(1,nTrial);
        sessionCostBin=ones(1,nTrial);
        % make bin
        quantileBenefit=quantile(sessionBenefitRating, (1/option.nBin):(1/option.nBin):((option.nBin-1)/option.nBin));
        quantileCost=quantile(sessionCostRating, (1/option.nBin):(1/option.nBin):((option.nBin-1)/option.nBin));
        for iQuantile=1:length(quantileBenefit)
             sessionBenefitBin(sessionBenefitRating>quantileBenefit(iQuantile))=sessionBenefitBin(sessionBenefitRating>quantileBenefit(iQuantile))+1;
             sessionCostBin(sessionCostRating>quantileCost(iQuantile))=sessionCostBin(sessionCostRating>quantileCost(iQuantile))+1;
        end
        benefitBin=[benefitBin sessionBenefitBin];
        costBin=[costBin sessionCostBin];
        % do linear fit
        regressor=glmfit([zscore(sessionBenefitRating)' zscore(sessionCostRating)'],sessionIsGoChoice', 'binomial', 'link', 'logit');
        subResult.firstLevelStat.betaBenefit(1,1,iSession)=regressor(2);
        subResult.firstLevelStat.betaCost(1,1,iSession)=regressor(3);
        regressor=glmfit([zscore(sessionBenefitRating)'],sessionIsGoChoice', 'binomial', 'link', 'logit');
        subResult.firstLevelStat.betaBenefitFromSeparateGLM(1,1,iSession)=regressor(2);
        regressor=glmfit([zscore(sessionCostRating)'],sessionIsGoChoice', 'binomial', 'link', 'logit');
        subResult.firstLevelStat.betaCostFromSeparateGLM(1,1,iSession)=regressor(2);
    end
    % fill subresult structure
    subResult.singleStat.binOfBenefitRatingValue = tools.tapply(benefitRating, {constant, benefitBin, subData.condition.sessionNumber}, @mean);
    subResult.singleStat.binOfBenefitChoice = tools.tapply(isGoChoice, {constant, benefitBin, subData.condition.sessionNumber}, @mean);
    subResult.singleStat.binOfCostRatingValue = tools.tapply(costRating, {constant, costBin, subData.condition.sessionNumber}, @mean);
    subResult.singleStat.binOfCostChoice = tools.tapply(isGoChoice, {constant, costBin, subData.condition.sessionNumber}, @mean);
    subResult.firstLevelStat.meanGo = tools.tapply(isGoChoice, {constant,constant,subData.condition.sessionNumber}, @mean);
    % fill result structure
    result.(submanipFieldName{iSubManip})= subResult;
end



