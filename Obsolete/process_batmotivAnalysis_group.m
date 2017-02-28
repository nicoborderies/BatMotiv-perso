function [design, statistics, comparison, correlations, dimensionality]= process_batmotivAnalysis_group( design )
%% define groups
group1 = 'CONTROL';
group2 = 'CITALOPRAM';

index = find(~ismember(design.Properties.VariableNames,{'subject_battery','subject','group','subject_neuropsyTable'})...
            & table2array(varfun(@isnumeric,design)));

%% descriptive statistics

statistics = table( nanmean(table2array(design(design.group==group1  | design.group==group2 ,index)),1)',...
                    nanstd(table2array(design(design.group==group1  | design.group==group2 ,index)),1)');
                
statistics.Properties.VariableNames = {['mean_' group1 '_' group2 ],['std_' group1 '_' group2]};
statistics.Properties.RowNames = design.Properties.VariableNames(index) ;    


%% ttest
        
        
[h,p,ci,stat] = ttest2(table2array(design(design.group==group1,index)),...
               table2array(design(design.group==group2,index)));
           
comparison = table( nanmean(table2array(design(design.group==group1  | design.group==group2 ,index)),1)',...
                    tools.sem(table2array(design(design.group==group1  | design.group==group2 ,index)),1)',...
                    nanmean(table2array(design(design.group==group2,index)),1)',...
                    tools.sem(table2array(design(design.group==group2,index)),1)',...
                    stat.tstat',...
                    p',...
                    h');
comparison.Properties.VariableNames = {['mean_' group1 ],['sem_' group1 ],['mean_' group2 ],['sem_' group2 ],'t_stat','pValue','H1'};
comparison.Properties.RowNames = design.Properties.VariableNames(index) ;    

                

%% correlation matrix
[correlations.population.rho,correlations.population.pValue] = corr(table2array(design(:,index)),'row','pairwise','type','Pearson');
[correlations.(group1).rho,correlations.(group1).pValue]  = corr(table2array(design(design.group==group1,index)),'row','pairwise','type','Pearson');
[correlations.(group2).rho,correlations.(group2).pValue]  = corr(table2array(design(design.group==group2,index)),'row','pairwise','type','Pearson');
sampleNames = {'population',group1,group2};
for iS = 1:numel(sampleNames)
    correlations.(sampleNames{iS}).rho = array2table(correlations.(sampleNames{iS}).rho);
    correlations.(sampleNames{iS}).rho.Properties.VariableNames = design.Properties.VariableNames(index) ;
    correlations.(sampleNames{iS}).rho.Properties.RowNames = design.Properties.VariableNames(index) ;
    correlations.(sampleNames{iS}).pValue = array2table(correlations.(sampleNames{iS}).pValue);
    correlations.(sampleNames{iS}).pValue.Properties.VariableNames = design.Properties.VariableNames(index) ;
    correlations.(sampleNames{iS}).pValue.Properties.RowNames = design.Properties.VariableNames(index) ;
end

%%  dimensionality reduction
% [dimensionality.coefficients,dimensionality.score,dimensionality.Variance,~,dimensionality.explainedVariance] = pca(table2array(design(:,index)),'Algorithm','eig','Rows','pairwise');
% 
% nComponent = 3;
% dimensionality.component = cell2table(cell(6,numel(dimensionality.coefficients(:,1))));
% % dimensionality.component.Properties.VariableNames = statistics.Properties.RowNames;
% dimensionality.component.Properties.RowNames = {'var component 1','var component 2','var component 3','x','xx','xxx'};
% for iComp = 1:nComponent
%     [coeff, ind ] = sort(dimensionality.coefficients(:,iComp));
%     coeffNames = cell2table(statistics.Properties.RowNames(ind)'); coeffNames.Properties.RowNames = {dimensionality.component.Properties.RowNames{iComp}};
%     dimensionality.component(iComp,:) = coeffNames;
% %     dimensionality.component(iComp+1,:) =coeff';
% end
% %% display
% f = figure;f.Color = 'w';
% scatter3(dimensionality.score(:,1),dimensionality.score(:,2),dimensionality.score(:,3),'filled');
% xlabel('1rst component');
% ylabel('2nd component');
% zlabel('3rd component');



end

