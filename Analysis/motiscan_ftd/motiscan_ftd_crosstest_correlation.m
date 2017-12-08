%% motiscan_ftd_crosstest_correlation


motiscanMetrics = [ stat.rating.mean_r , stat.rating.mean_e,...
                    stat.choice.choice_accuracy_r , stat.choice.choice_accuracy_e , ...
                    stat.weight.acceptRate ,...
                    stat.discount.patientChoice,...
                    stat.grip.mean_nfpeak,stat.grip.incentive_fpeak,...
                    stat.learning.correct_gain , stat.learning.correct_loss ]; 

group = stat.rating.group;

% correlation matrix
corrmat = cell(1,2);
for ig = 1:2
   
    corrmat{ig} = corr(motiscanMetrics(group==ig,:),'type','Spearman','rows','pairwise');
    
end

%%
% display
metricLabels = {'rating_r','rating_e',...
                'choice_accuracy_r','choice_accuracy_e',...
                'acceptRate_re',...
                'patientChoice_rd',...
                'mean_fpeak','incentive_fpeak',...
                'correct_gain','correct_loss'};
groupLabel = {'CONTROL','FTD'};

for ig = 1:2
   
    f1 = figure;
    h = heatmap(corrmat{ig},[],[],1,'TextColor',[1 1 1]);
    colormap jet; caxis([-1 1]);colorbar;
    ax = gca;
    title([ groupLabel{ig} ' : metrics correlation matrix']);
    ax.XAxisLocation = 'bottom';
    ax.TickLength = [0 0];
    ax.XTick = [1:size(corrmat{ig},1)];
    ax.YTick = [1:size(corrmat{ig},1)];
    ax.XTickLabel = metricLabels;
    ax.YTickLabel = metricLabels;
    ax.XTickLabelRotation = 90;
    ax.TickLabelInterpreter = 'tex';
    box on;
    f1.Position = [ 100 100 900 900];
    set_all_properties('FontName','Arial Narrow','FontWeight','normal','FontSize',16,...
                        'LineWidth',1.5,'Interpreter','tex');
                    
end