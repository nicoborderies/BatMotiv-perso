%% hiddenvalue_model_identifiability
%
% nicolas borderies - 09/2017

clc;
clear all;
close all;

%% initialize

% Potential parameters
paramRange.alpha = [0.1 1 5];
paramRange.theta = [0.1 1 5];
paramRange.mu  = [0:1:5];
paramRange.sigma = [0.1 3 5];

% Potential observation combinations
get_choice = [0 1];
get_accept = [0 1];
get_responsetime = [0 1];
precisionRange = 10.^([0:4]);
  

% Potential vectors of parameters/outputs combinations
PARAM = combvec(paramRange.alpha,paramRange.theta,paramRange.mu,paramRange.sigma)';
OUTPUTS = combvec(get_choice,get_accept,get_responsetime,precisionRange)';

% set parameters vectors
MU = 1;
SIGMA = 1;
ALPHA = 1;
THETA = 1;
% MU = PARAM(:,1);
% SIGMA = PARAM(:,2);
% ALPHA = PARAM(:,3);
% THETA = PARAM(:,4);

% set ouptus combinations
CHOICE = OUTPUTS(:,1);
ACCEPT = OUTPUTS(:,2);
RT = OUTPUTS(:,3);
PRECISION = OUTPUTS(:,4);

% set number of simulations
% nsim = size(PARAM,1);
% nsim = size(OUTPUTS,1);
nsim = 30;
ncombinations = size(OUTPUTS,1);


%% simulations
result = cell(1,nsim);
estimates = cell(1,nsim);

for isim = 1:nsim
    
    clc;
    fprintf(' simulation: %d/%d \n',isim,nsim);

    % parameter-space exploration
%     param.mu = MU(isim);
%     param.sigma =  SIGMA(isim);
%     param.alpha =  ALPHA(isim);
%     param.theta = THETA(isim);
    
    % parameter radomization
    param.mu = rand(1)*max(paramRange.mu) + min(paramRange.mu) ;
    param.sigma =  rand(1)*max(paramRange.sigma) + min(paramRange.sigma) ;
    param.alpha =  rand(1)*max(paramRange.alpha) + min(paramRange.alpha) ;
    param.theta = rand(1)*max(paramRange.theta) + min(paramRange.theta) ;
    
    for icombinations = 1:ncombinations
        
        fprintf(' combination: %d/%d \n',icombinations,ncombinations);

        
        % options
        param.get_choice = CHOICE(icombinations);
        param.get_accept = ACCEPT(icombinations);
        param.get_responsetime = RT(icombinations);
        param.precision = PRECISION(icombinations);

        % index
        iloop = icombinations + (isim-1)*ncombinations;
        
        % generative predictions
        [ result{iloop} ] = predict_hiddenvalue( param );

        % inferential estimation
        [estimates{iloop}] = estimate_hiddenvalue(result{iloop}.inputs,result{iloop}.outputs,param);
    end

end

% store results
save('hiddenvalue_model_identifiability_informationaldependance.mat',...
     'PARAM','OUTPUTS','result','estimates');

%% observational dependance analysis

% - extract estimates
ncondition = 4;
nparam = 4;
ncombination = 40;
nsample = numel(result);
generativeParam = nan(nsample,nparam);
estimatedParam = nan(nsample,nparam);
observationConditions = nan(nsample,ncondition);

for isample=1:nsample
    generativeParam(isample,:) = [ result{isample}.param.alpha , result{isample}.param.theta , result{isample}.param.mu , result{isample}.param.sigma ];
    estimatedParam(isample,:) = [ estimates{isample}.param.alpha , estimates{isample}.param.theta , estimates{isample}.param.mu , estimates{isample}.param.sigma ];
    observationConditions(isample,:) = [ result{isample}.param.get_choice , result{isample}.param.get_accept , result{isample}.param.get_responsetime , result{isample}.param.precision];
end

% - Analyse across observation conditions
observationCode = findgroups(array2table(observationConditions));
observationSet = unique(observationCode)';
identifiability = cell(1,numel(observationSet));

for cond = observationSet
    
    % - define
    identifiability{cond} = struct;
    % - selection of a subset
    subset = (observationCode==cond);
    % - compute confusion matrix - full correlation
    % confusionmat = corr(estimatedParam,generativeParam);
    % - compute confusion matrix - partial correlation
    confusionmat = nan(nparam);
    for i=1:nparam
        for j=1:nparam
    %         confusionmat(i,j) = partialcorr(estimatedParam(:,i),generativeParam(:,j),estimatedParam(:,setdiff([1:nparam],i)));
            beta = glmfit(zscore([ estimatedParam(subset,i) , estimatedParam(subset,setdiff([1:nparam],i)) ]),zscore(generativeParam(subset,j)),'normal','constant',0);
            confusionmat(i,j) = beta(1);
        end
    end
    determinationmat = confusionmat.^2;
    % - accuracy
    accuracy = nan(1,nparam);
    for i=1:nparam
        accuracy(i) = determinationmat(i,i);
%         accuracy(i) = mean([determinationmat(i,i),mean(1-[determinationmat(setdiff([1:nparam],i),i)])]);
    end
    % - store
    identifiability{cond} = struct('confusionmat',confusionmat,...
                                   'determinationmat',determinationmat,...
                                   'accuracy',accuracy);
end

% store results
save('hiddenvalue_model_identifiability_informationaldependance.mat',...
     'identifiability',...
     '-append');

%% display
paramLabels = {'\alpha','\theta','\mu','\sigma'};

% parameters confusion matrix
confusionmat = identifiability{end}.confusionmat;
f1 = figure;
h = heatmap(confusionmat,[],[],1,'TextColor',[1 1 1]);
% h = heatmap(determinationmat,[],[],1,'TextColor',[1 1 1]);
colormap jet; caxis([-1 1]);colorbar;
ax = gca;
title('parameters confusion matrix');
ax.XLabel.String = 'generative parameters';
ax.YLabel.String = 'estimated parameters';
ax.XAxisLocation = 'top';
ax.TickLength = [0 0];
ax.XTick = [1:4];
ax.YTick = [1:4];
ax.XTickLabel = paramLabels;
ax.YTickLabel = paramLabels;
ax.TickLabelInterpreter = 'tex';
box on;
f1.Position = [ 100 100 600 600];
set_all_properties('FontName','Arial Narrow','FontWeight','normal','FontSize',16,...
                    'LineWidth',1.5,'Interpreter','tex');
                
%% accuracy display

f2 = figure;
for iparam = 1:nparam
    subplot(nparam,1,iparam); hold on;
    for cond = observationSet
        bar(cond,identifiability{cond}.accuracy(iparam),'FaceColor',[1 1 1]*0.5);
    end
    ax = gca;
    ax.TickLength = [0 0];
    ax.XTick = [];
    ylim([0 1]);
    xlim([0 numel(observationSet)+1]);
    box off;
    title(paramLabels{iparam},'FontSize',16);
    ylabel(' estimation accuracy (%)','FontSize',14);
    if iparam==nparam
        ax.XLabel.String = 'observational conditions';
        ax.XLabel.FontSize = 16;
    end
end
f2.Position = [ 100 100 600 900];
set_all_properties('FontName','Arial Narrow','FontWeight','normal',...
                    'LineWidth',1.5,'Interpreter','tex');
                
%%
conditionLabels = {'choice_2O','choice_1O','rt','precision'};
conditionTab = array2table(OUTPUTS','RowNames',conditionLabels);
f3 = figure;
heatmap(conditionTab{:,:},[],[],1,'TextColor',[1 1 1]*0.5);
cmap = colormap('gray'); colormap(1-cmap); caxis([0 1]);
ax = gca;
ax.YTick = [1:numel(conditionLabels)];
ax.YTickLabel = conditionLabels;
box on;
f3.Position = [ 100 100 600 100];
set_all_properties('FontName','Arial Narrow','FontWeight','normal','FontSize',16,...
                    'LineWidth',1.5,'Interpreter','tex');


