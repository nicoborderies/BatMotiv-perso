function [result] = process_data_rating_modelfree(dataStructure, option)
% process_data_rating_modelfree do the model free analysis for the task rating of the MBB Battery for a single subject
% [data] = process_data_rating_modelfree(dataStructure)
%
% IN
%       - data structure from the load_data_rating function
%       - No option here
% OUT
%       - results containing one field by submanip. Each submanip contains:
%       -------------------------------------------------------------------------------------------------------------
%           * singleStat : statistics usefull for a single subject (e.g. to plot) but won't be passed to the second level... ex variance, distribution etc
%               - percentileRating     : [.025 .25 .50 .75 .975] percentiles of ratings (1*5*nSession)
%           * firstLevelStat : statistics easy to pass to the second level (e.g for t.tests)
%               - meanRating                   : mean value of ratings (1*1*nSession)
%               - meanRatingByItemSubtype      : mean value of ratings by item Subtype (2*1*nSession)
%               - varianceRating               : variance value of ratings (1*1*nSession)
%               - skewnessRating               : skewness value of ratings (assymetry metric) (1*1*nSession)
%               - kurtosisRating               : kurtosis value of ratings (peakedness metric)(1*1*nSession)


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
    subData=dataStructure.(submanipFieldName{iSubManip});
    constant=ones(size(subData.condition.sessionNumber));
    for iSession = unique(subData.condition.sessionNumber)
        subResult.singleStat.percentileRating(1,:,iSession) = quantile(subData.behavior.rating(subData.condition.sessionNumber==iSession),  [.025 .25 .50 .75 .975]);
    end
    subResult.firstLevelStat.meanRating = tools.tapply(subData.behavior.rating, {constant,constant,subData.condition.sessionNumber}, @mean);
    subResult.firstLevelStat.varianceRating = tools.tapply(subData.behavior.rating, {constant,constant,subData.condition.sessionNumber}, @var);
    subResult.firstLevelStat.skewnessRating = tools.tapply(subData.behavior.rating, {constant,constant,subData.condition.sessionNumber}, @skewness);
    subResult.firstLevelStat.kurtosisRating = tools.tapply(subData.behavior.rating, {constant,constant,subData.condition.sessionNumber}, @kurtosis);
    subResult.firstLevelStat.meanRatingByItemSubtype = tools.tapply(subData.behavior.rating, {constant,subData.condition.itemSubtype,subData.condition.sessionNumber}, @mean);
    
    result.(submanipFieldName{iSubManip})= subResult;
end


