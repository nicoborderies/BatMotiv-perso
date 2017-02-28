function [design, statistics, correlations , dimensionality]= process_batmotiv_group( design, data, misc )
%% data selection
groupList = {'CONTROL','CITALOPRAM'};

% first selection 
    % select subject sessions
        factor = {'treatment'};iFactor = 2;
%         indexRow = find(ismember(design.group,groupList)...
%                     & design.(factor{1}) ~= iFactor...
%                     & design.sessionNumber ~= 2);
                
        indexRow = (misc.treatment~=2);

    % select variables
        % select numeric variables
        indexVar = find(table2array(varfun(@isnumeric,design))...
                    & ~ismember(design.Properties.VariableNames,{'subject','sessionNumber','treatment'}) );

        % select scalar variables
        indexDoubleVar = find(table2array(varfun(@size,design))'==2)/2;
        indexVar = indexVar(~ismember(indexVar,indexDoubleVar));

        design = design(indexRow,indexVar);
        misc   = misc(misc.(factor{1}) ~= iFactor & misc.sessionNumber ~= 2,:);


% second selection

    indexVar = 1:width(design);
    
    % select defined variables
    indexUndef = find(sum(table2array(varfun(@isnan,design))==0)==0);
    indexVar = indexVar(~ismember(indexVar,indexUndef));

    % select free (not-fixed) variables
    fixedVar = {'learning_temperature'};
    indexFixed = find( ismember(design.Properties.VariableNames,fixedVar));
    indexVar = indexVar(~ismember(indexVar,indexFixed));
    
    design = design(:,indexVar);


%% descriptive statistics

statistics = table(table2array(varfun(@nanmean,design))' , table2array(varfun(@nanstd,design))');
statistics.Properties.VariableNames = {['mean'],['std']};
statistics.Properties.RowNames = design.Properties.VariableNames ;    

misc = data.battery.misc ;
mat  = table2array(misc(:,3:end));
sum(mat~=0,1);
misc{height(misc)+1,2:end} = [  height(misc) , sum(mat~=0,1) ];
misc{end,1} = nominal(numel(unique(misc.group(~isundefined(misc.group))))) ;

% export
writetable(statistics,'descriptiveStatistics_batmotiv_population.xlsx','WriteRowNames',1,'Sheet',1);
writetable(misc,'descriptiveStatistics_batmotiv_population.xlsx','WriteRowNames',1,'Sheet',2);



%% ttest
        
%         
% [h,p,ci,stat] = ttest2(table2array(design(design.group==group1,index)),...
%                table2array(design(design.group==group2,index)));
%            
% comparison = table( nanmean(table2array(design(design.group==group1  | design.group==group2 ,index)),1)',...
%                     tools.sem(table2array(design(design.group==group1  | design.group==group2 ,index)),1)',...
%                     nanmean(table2array(design(design.group==group2,index)),1)',...
%                     tools.sem(table2array(design(design.group==group2,index)),1)',...
%                     stat.tstat',...
%                     p',...
%                     h');
% comparison.Properties.VariableNames = {['mean_' group1 ],['sem_' group1 ],['mean_' group2 ],['sem_' group2 ],'t_stat','pValue','H1'};
% comparison.Properties.RowNames = design.Properties.VariableNames(index) ;    

                

%% correlation matrix

[correlations.rho,correlations.pValue] = corr(table2array(design),'row','pairwise','type','Pearson');
correlations.rho = array2table(correlations.rho);
correlations.rho.Properties.VariableNames = design.Properties.VariableNames;
correlations.rho.Properties.RowNames      = design.Properties.VariableNames ;
correlations.pValue = array2table(correlations.pValue);
correlations.pValue.Properties.VariableNames = design.Properties.VariableNames ;
correlations.pValue.Properties.RowNames = design.Properties.VariableNames;

% ordonnate correlations
pV = table2array(correlations.pValue); 
pV2 = nan(1,numel(pV));
pV2(find(pV)) = pV(find(pV));
[pV3,iP] = sort(pV2);
pV4 = pV3(pV3<=0.05);
iP4 = iP(pV3<=0.05);

orderCorr = table( categorical(ones(numel(pV4),1)) , categorical(ones(numel(pV4),1)) ,nan(numel(pV4),1) ,  nan(numel(pV4),1) );
orderCorr.Properties.VariableNames = {'var1','var2','rho','pValue'};
for iCor = 1:numel(pV4)
    [l,c] = find(correlations.pValue{:,:}==pV4(iCor));
    orderCorr.var1(iCor) = correlations.pValue.Properties.VariableNames{l(1)};
    orderCorr.var2(iCor) = correlations.pValue.Properties.VariableNames{c(1)};
    orderCorr.rho(iCor) = correlations.rho{l(1),c(1)};
    orderCorr.pValue(iCor) = correlations.pValue{l(1),c(1)};
end
% %  display
% displayCorrelations(correlations,design)


%%  dimensionality reduction

[dimensionality.coefficients,dimensionality.score,dimensionality.Variance,~,dimensionality.explainedVariance] = pca(table2array(design),'Algorithm','eig','Rows','pairwise');
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

% % display
% f = figure;f.Color = 'w';
% scatter3(dimensionality.score(:,1),dimensionality.score(:,2),dimensionality.score(:,3),'filled');
% xlabel('1rst component');
% ylabel('2nd component');
% zlabel('3rd component');

%% model comparison

% data = 'choice';
% data = 'weight';
data = 'battery';
subSelect = (misc.treatment~=2);
% subSelect = (isnan(misc.treatment)); % control group
% subSelect = (isnan(misc.sessionNumber) & misc.treatment~=2); % citalopram group
misc2 = misc(subSelect,:);

% options.families = {[1:2:23],[2:2:24]};
% options.families = {[1:24]};


logE = [];
analysis = unique(misc2.analysisNumber);
for iA = 1:numel(analysis)
    logE = [logE ; misc2.(['logE_' data ])(misc2.analysisNumber==iA)'];
end
logE = logE(:,~isnan(logE(1,:)));


[posterior,out] = VBA_groupBMC(logE);
% [posterior,out] = VBA_groupBMC(logE,options);


end

