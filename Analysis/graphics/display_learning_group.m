function [fig] = display_learning_2(dataTable,descriptiveTable,inferentialTable,option)

% plot_learning plot the stochastic reinforcement learning task of the MBB Battery at
% the group level, under placebo vs. citalopram

if nargin <3;
    fig = [];end



%% FIGURES %%
%__________________________________________________________________________
% Figure Specifications
colorTreatment =[0, 0 ,0 ; 1, 0, 1];
colorSession = [0, 0 ,0 ; 0, 0, 1];
colorGroup = {'b','r'};
title_size = 12;
axis_size=14;tick_size=10;
legend_size=14;
xSize= 40; ySize= 20;
x_text = 80; y_text = 80; y_shift = 10; text_color = '{0 0 0}';
star_size = 14;marker_size=0.3;




%----------------------------------------- Figure 1 ----------------------------------------------------%
% line plots / bar plots: probability of optimal choice as a function of session evolution / valence %
perfPlotLearning = figure;set(perfPlotLearning, 'color', 'white');
perfPlotLearning.Units = 'centimeters'; perfPlotLearning.Position = [ 1 1 30 14];


% Performance Evolution
subplot(121);
hold on;box on;
    
    % data
    for iGroup = 1:2
        eval([ 'data' num2str(iGroup) ' = groupData.(groupNames{groupCompare(iGroup)}).group.learning.gainLearning ;' ]);
        eval([ 'errorData' num2str(iGroup) ' = groupData.(groupNames{groupCompare(iGroup)}).group.learning.errorGainLearning ;' ]);
        eval([ 'prediction' num2str(iGroup) ' = groupData.(groupNames{groupCompare(iGroup)}).group.learning.predictedGainLearning ;' ]);
    end

    % random/matching levels
    for proba = [0.5 0.66]
        plot([0:numel(data1)+1],ones(numel(data1)+2,1)*proba,'Color',colorTreatment(1,:),'Color',[0.6 0.6 0.6],'LineWidth',1,'LineStyle','--');
    end
    
    % evolution for gain
    [h,hp] = boundedline(0:numel(data1), [0.5 data1], [0 errorData1], 'alpha');
    set(h,'Color',colorGroup{1},'LineWidth',2);set(hp,'FaceColor',colorGroup{1});
    h1=scatter(1:numel(data1), data1,80,colorGroup{1},'filled');
%     h = errorbar(0:numel(data1), [0.5 data1], [0 errorData1]);
%     set(h,'Color','b','LineWidth',2); h.LineStyle='none';
    plot(0:numel(data1),[0.5 prediction1],'b--');
    
    [h,hp] = boundedline(0:numel(data2), [0.5 data2], [0 errorData2], 'alpha');
    set(h,'Color',colorGroup{2},'LineWidth',2);set(hp,'FaceColor',colorGroup{2});
    h2=scatter(1:numel(data2), data2,80,colorGroup{2},'filled');
%     h = errorbar(0:numel(data1), [0.5 data1], [0 errorData1]);
%     set(h,'Color','b','LineWidth',2); h.LineStyle='none';
    plot(0:numel(data2),[0.5 prediction2],'r--');

    % axis properties
        ax = get(gca);
    set(gca,'YLim',[0 1],...
    'XLim',[0 6],...
    'XTick',[0 1 2 3 4 5 6],...
    'TickLength',[0 0],...
    'YTick',[0 0.25 0.5 0.75 1],...
    'YTickLabel',{'0';'25';'50';'75';'100'});
    ax.XTickLabel = {'0';'4';'8';'12';'16';'20';'24'};

    xlabel('numéro d''essai','FontSize',14);
    ylabel({'choix corrects (%)'},'FontSize',14);
    title('apprentissage aux gains','FontSize',14);
    ylim([0 1]);
    
     % legend
     legend([h1 h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
    

subplot(122);
hold on;box on;
    
    % data
    for iGroup = 1:2
        eval([ 'data' num2str(iGroup) ' = groupData.(groupNames{groupCompare(iGroup)}).group.learning.lossLearning ;' ]);
        eval([ 'errorData' num2str(iGroup) ' = groupData.(groupNames{groupCompare(iGroup)}).group.learning.errorLossLearning ;' ]);
        eval([ 'prediction' num2str(iGroup) ' = groupData.(groupNames{groupCompare(iGroup)}).group.learning.predictedLossLearning ;' ]);
    end

    % random/matching levels
    for proba = [0.5 0.83]
        plot([0:numel(data1)+1],ones(numel(data1)+2,1)*proba,'Color',colorTreatment(1,:),'Color',[0.6 0.6 0.6],'LineWidth',1,'LineStyle','--');
    end
    
    % evolution for loss
    [h,hp] = boundedline(0:numel(data1), [0.5 data1], [0 errorData1], 'alpha');
    set(h,'Color',colorGroup{1},'LineWidth',2);set(hp,'FaceColor',colorGroup{1});
    h1=scatter(1:numel(data1), data1,80,colorGroup{1},'filled');
%     h = errorbar(0:numel(data1), [0.5 data1], [0 errorData1]);
%     set(h,'Color','b','LineWidth',2); h.LineStyle='none';
    plot(0:numel(data1),[0.5 prediction1],'b--');
    
    [h,hp] = boundedline(0:numel(data2), [0.5 data2], [0 errorData2], 'alpha');
    set(h,'Color',colorGroup{2},'LineWidth',2);set(hp,'FaceColor',colorGroup{2});
    h2=scatter(1:numel(data2), data2,80,colorGroup{2},'filled');
%     h = errorbar(0:numel(data1), [0.5 data1], [0 errorData1]);
%     set(h,'Color','b','LineWidth',2); h.LineStyle='none';
    plot(0:numel(data2),[0.5 prediction2],'r--');

    % axis properties
        ax = get(gca);
    set(gca,'YLim',[0 1],...
    'XLim',[0 6],...
    'XTick',[0 1 2 3 4 5 6],...
    'TickLength',[0 0],...
    'YTick',[0 0.25 0.5 0.75 1],...
    'YTickLabel',{'0';'25';'50';'75';'100'});
    ax.XTickLabel = {'0';'4';'8';'12';'16';'20';'24'};

    xlabel('numéro d''essai','FontSize',14);
    ylabel({'choix corrects (%)'},'FontSize',14);
    title('apprentissage aux pertes','FontSize',14);
    ylim([0 1]);
    
     % legend
     legend([h1 h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
    
% store figure handles
fig = perfPlotLearning;

end