function [fig] = display_learning_2(dataTable,descriptiveTable,inferentialTable,option)
%% Specifications
    

%% Figure parameters

colorTreatment =[0, 0 ,0 ; 1, 0, 1];
colorSession = [0, 0 ,0 ; 0, 0, 1];
colorGroup = {'b','r'};
title_size = 12;
axis_size=14;tick_size=10;
legend_size=14;
xSize= 40; ySize= 20;
x_text = 80; y_text = 80; y_shift = 10; text_color = '{0 0 0}';
star_size = 14;marker_size=0.3;

nbin=24;



%% Fig 1: choice dependency to benefits, costs, subjective value
%--------------------------------
% define fig handle

learningFig1 = figure;set(learningFig1, 'color', 'white');
learningFig1.Units = 'centimeters'; learningFig1.Position = [ 1 1 30 15];


 % plot selection
        % -------------- errorbar of performance evolution
            subplot(1,1,1);  hold on;
            % data to plot
                KR = nanmean(inferentialTable.(['kR_learning']));
                stdKR = nanstd(inferentialTable.(['kR_learning']));
                KP = nanmean(inferentialTable.(['kP_learning']));
                stdKP = nanstd(inferentialTable.(['kP_learning']));
                ALPHA = nanmean(inferentialTable.(['alpha']));
                stdALPHA= nanstd(inferentialTable.(['alpha']));
            
                for iValence = [-1,1]
                    selection =  ismember(dataTable.group,{'CONTROL','CITALOPRAM'})...
                                & dataTable.pairValence~=0 ...
                                & dataTable.pairValence==iValence; % select group
                    subjects = dataTable.subjectNumber; subList = unique(subjects);
                    TRIALS = []; CHOICE=[];PREDICTION = [];
                    for iSub = 1:numel(subList)
                         trials = dataTable.trialNumberByValence(selection & subjects==subList(iSub));   
                         choice = dataTable.isOptimalChoice(selection & subjects==subList(iSub));   
    %                      prediction = dataTable.predicted_goChoice(selection & subjects==subList(iSub));  
                         TRIALS(iSub,:) = tools.tapply(trials,{trials},@nanmean,{'continuous'},[nbin]);
                         CHOICE(iSub,:) = tools.tapply(choice,{trials},@nanmean,{'continuous'},[nbin]);
%                          PREDICTION(iSub,:) = tools.tapply(prediction,{trials},@nanmean,{'continuous'},[nbin]);
                    end

                
                % plot
                    switch iValence
                        case 1
                            [h,hp] = boundedline( nanmean(TRIALS,1) , nanmean(CHOICE,1) , sem(CHOICE,1) , 'alpha' );   
                            set(h,'Color','g','LineWidth',2);set(hp,'FaceColor','g');
                            h.LineStyle = 'none'; 
                            scatter( nanmean(TRIALS,1) , nanmean(CHOICE,1) , 'filled','g' );
                        case -1
                            [h,hp] = boundedline( nanmean(TRIALS,1) , nanmean(1-CHOICE,1) , sem(CHOICE,1) , 'alpha' );   
                            set(h,'Color','r','LineWidth',2);set(hp,'FaceColor','r');
                            h.LineStyle = 'none'; 
                            scatter( nanmean(TRIALS,1) , nanmean(1-CHOICE,1) , 'filled','r' );
                    end

                end
                
            % legend
                xx = xlim; xx = 1; yy = 1;
                xlabel(texlabel([ ' trial number' ]),'FontSize',14);
                ylabel(' correct choice (%)','FontSize',14);
                text(xx,yy,texlabel(['k_R = ' num2str(round(KR,2)) ' (+/- ' num2str(round(stdKR,2)) ')']));
                text(xx,yy-0.1,texlabel(['k_P = ' num2str(round(KP,2)) ' (+/- ' num2str(round(stdKP,2)) ')']));
                text(xx,yy-0.2,texlabel(['alpha = ' num2str(round(ALPHA,2)) ' (+/- ' num2str(round(stdALPHA,2)) ')']));


    %             legend([h h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
            % axis 
                ax = get(gca);box off; ax.FontName = 'FixedWidth';ax.LabelFontSizeMultiplier = 0.9;
                set(gca,...
                    'XTick',[1:nbin],...
                    'YLim',[0 1],...
                    'TickLength',[0 0]);

% subplot(122);
% hold on;box on;
%     
%     % data
%     for iGroup = 1:2
%         eval([ 'data' num2str(iGroup) ' = groupData.(groupNames{groupCompare(iGroup)}).group.learning.lossLearning ;' ]);
%         eval([ 'errorData' num2str(iGroup) ' = groupData.(groupNames{groupCompare(iGroup)}).group.learning.errorLossLearning ;' ]);
%         eval([ 'prediction' num2str(iGroup) ' = groupData.(groupNames{groupCompare(iGroup)}).group.learning.predictedLossLearning ;' ]);
%     end
% 
%     % random/matching levels
%     for proba = [0.5 0.83]
%         plot([0:numel(data1)+1],ones(numel(data1)+2,1)*proba,'Color',colorTreatment(1,:),'Color',[0.6 0.6 0.6],'LineWidth',1,'LineStyle','--');
%     end
%     
%     % evolution for loss
%     [h,hp] = boundedline(0:numel(data1), [0.5 data1], [0 errorData1], 'alpha');
%     set(h,'Color',colorGroup{1},'LineWidth',2);set(hp,'FaceColor',colorGroup{1});
%     h1=scatter(1:numel(data1), data1,80,colorGroup{1},'filled');
% %     h = errorbar(0:numel(data1), [0.5 data1], [0 errorData1]);
% %     set(h,'Color','b','LineWidth',2); h.LineStyle='none';
%     plot(0:numel(data1),[0.5 prediction1],'b--');
%     
%     [h,hp] = boundedline(0:numel(data2), [0.5 data2], [0 errorData2], 'alpha');
%     set(h,'Color',colorGroup{2},'LineWidth',2);set(hp,'FaceColor',colorGroup{2});
%     h2=scatter(1:numel(data2), data2,80,colorGroup{2},'filled');
% %     h = errorbar(0:numel(data1), [0.5 data1], [0 errorData1]);
% %     set(h,'Color','b','LineWidth',2); h.LineStyle='none';
%     plot(0:numel(data2),[0.5 prediction2],'r--');
% 
%     % axis properties
%         ax = get(gca);
%     set(gca,'YLim',[0 1],...
%     'XLim',[0 6],...
%     'XTick',[0 1 2 3 4 5 6],...
%     'TickLength',[0 0],...
%     'YTick',[0 0.25 0.5 0.75 1],...
%     'YTickLabel',{'0';'25';'50';'75';'100'});
%     ax.XTickLabel = {'0';'4';'8';'12';'16';'20';'24'};
% 
%     xlabel('numéro d''essai','FontSize',14);
%     ylabel({'choix corrects (%)'},'FontSize',14);
%     title('apprentissage aux pertes','FontSize',14);
%     ylim([0 1]);
%     
%      % legend
%      legend([h1 h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
%     
%% Save fig handle
    fig{1} = weightFig1;  

end