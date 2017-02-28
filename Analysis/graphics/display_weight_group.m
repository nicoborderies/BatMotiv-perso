function [fig] = display_weight_group(groupNames,groupCompare,groupResult,groupData)
%% Specifications
    subTask = {'writtenReward','writtenPunishment','writtenEffort'};
    subTaskName = {'reward (text)','punishment (text)','effort (text)'};
    subTaskAcronym = {'R','P','E'};
    subTaskSubCategories = {'alimentary','nonAlimentary';
                            'sensory','nonSensory';
                            'motor','cognitive'};


%% Figure parameters
    colorTreatment =[0, 0 ,0 ; 1, 0, 1];
    colorSession = [0, 0 ,0 ; 0, 0, 1];
    title_size = 12;
    axis_size=14;tick_size=10;
    legend_size=14;
    xSize= 40; ySize= 20;
    x_text = 80; y_text = 80; y_shift = 10; text_color = '{0 0 0}';
    star_size = 14;marker_size=0.3;
    
    nbin=10;

%% Fig 1: choice dependency to deltaV
%--------------------------------
% define fig handle
    choiceFig1 = figure;set(choiceFig1, 'color', 'white');
    choiceFig1.Units = 'centimeters'; choiceFig1.Position = [ 1 1 35 10];

% iteration across subtask
    for iSubTask = 1:numel(subTask)
        
        % plot selection
        % -------------- errorbar of choice
            subplot(1,3,iSubTask);  hold on;
            % data to plot
                selection = ismember(dataTable.dimension,subTask{iSubTask})...       % select subtask
                            & ismember(dataTable.group,{'CONTROL','CITALOPRAM'}) ; % select group
                subjects = dataTable.subjectNumber; subList = unique(subjects);
                for iSub = 1:numel(subList)
                     dv = dataTable.differenceRatingValue(selection & subjects==subList(iSub));   
                     choice = dataTable.sideChoice(selection & subjects==subList(iSub));   choice = (choice+1)./2;
                     k = inferentialTable.(['k' subTaskAcronym{iSubTask}])(inferentialTable.subject==subList(iSub));
                     prediction = sig(k(1)*(dv/100));
                     DV(iSub,:) = tools.tapply(dv,{dv},@nanmean,{'continuous'},[nbin]);
                     CHOICE(iSub,:) = tools.tapply(choice,{dv},@nanmean,{'continuous'},[nbin]);
                     PREDICTION(iSub,:) = tools.tapply(prediction,{dv},@nanmean,{'continuous'},[nbin]);
                end
                K = nanmean(inferentialTable.(['k' subTaskAcronym{iSubTask}]));
                stdK = nanstd(inferentialTable.(['k' subTaskAcronym{iSubTask}]));
                ACCURACY = nanmean(inferentialTable.(['congruentChoice' subTaskAcronym{iSubTask}]));
                stdACCURACY = nanstd(inferentialTable.(['congruentChoice' subTaskAcronym{iSubTask}]));

                
            % plot
                h = errorbar( nanmean(DV,1) , nanmean(CHOICE,1) , sem(CHOICE,1) );                
                scatter( nanmean(DV,1) , nanmean(CHOICE,1) , 'filled','k' );
                plot( nanmean(DV,1) , nanmean(PREDICTION,1) , '-k' );
                h.LineStyle = 'none';  h.LineWidth = 1.5; h.Color = [0 0 0]; 
                
            % legend
                xx = xlim; xx = -95; yy = 0.9;
                xlabel(texlabel([ 'delta V (rating)' ]),'FontSize',14);
                if iSubTask==1; ylabel('choice (%)','FontSize',14);end
                title(subTaskName{iSubTask},'FontSize',14);
                text(xx,yy,texlabel(['k_' subTaskAcronym{iSubTask} ' = ' num2str(round(K,2)) ' (+/- ' num2str(round(stdK,2)) ')']));
                text(xx,yy-0.1,texlabel(['Accuracy_' subTaskAcronym{iSubTask} ' = ' num2str(round(ACCURACY,2)) ' (+/- ' num2str(round(stdACCURACY,2)) ')']));

    %             legend([h h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
            % axis 
                ax = get(gca);box off; ax.FontName = 'FixedWidth';ax.LabelFontSizeMultiplier = 0.9;
                set(gca,...
                    'XLim',[-100 100],...
                    'YLim',[0 1],...
                    'TickLength',[0 0]);
    end
    
%% Fig 2: RT dependency to deltaV
%--------------------------------
% define fig handle
    choiceFig2 = figure;set(choiceFig2, 'color', 'white');
    choiceFig2.Units = 'centimeters'; choiceFig2.Position = [ 1 1 35 10];

% iteration across subtask
    for iSubTask = 1:numel(subTask)
        
        % plot selection
        % -------------- errorbar of choice
            subplot(1,3,iSubTask);  hold on;
            % data to plot
                selection = ismember(dataTable.dimension,subTask{iSubTask})...       % select subtask
                            & ismember(dataTable.group,{'CONTROL','CITALOPRAM'}) ; % select group
                subjects = dataTable.subjectNumber; subList = unique(subjects);
                for iSub = 1:numel(subList)
                     dv = abs(dataTable.differenceRatingValue(selection & subjects==subList(iSub)));   
                     rt = dataTable.rt(selection & subjects==subList(iSub));   
                     DV(iSub,:) = tools.tapply(dv,{dv},@nanmean,{'continuous'},[nbin]);
                     RT(iSub,:) = tools.tapply(rt,{dv},@nanmean,{'continuous'},[nbin]);
                end

            % plot
                h = errorbar( nanmean(DV,1) , nanmean(RT,1) , sem(RT,1) );                
                scatter( nanmean(DV,1) , nanmean(RT,1) , 'filled','k' );
                h.LineStyle = '-';  h.LineWidth = 1.5; h.Color = [0 0 0]; 
                
            % legend
                xlabel(texlabel([ '|deltaV| (rating)' ]),'FontSize',14);
                if iSubTask==1; ylabel('response time (sec)','FontSize',14);end
                title(subTaskName{iSubTask},'FontSize',14);

    %             legend([h h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
            % axis 
                ax = get(gca);box off; ax.FontName = 'FixedWidth';ax.LabelFontSizeMultiplier = 0.9;
                set(gca,...
                    'TickLength',[0 0]);
    end
    
%% Save fig handle
    fig{1} = choiceFig1;    
    fig{2} = choiceFig2;    
    
    

% %% FIGURES %%
% %__________________________________________________________________________
% % Figure Specifications
% title_size = 12;
% axis_size=14;tick_size=10;
% legend_size=14;
% xSize= 40; ySize= 20;
% x_text = 80; y_text = 80; y_shift = 10; text_color = '{0 0 0}';
% star_size = 14;marker_size=0.3;
% map(1,:) = [1 0 0]; % colormap
% for nContrast = 2:32
%     map(nContrast,:) = [map(nContrast-1,1), map(nContrast-1,2)+(1/31), map(nContrast-1,3)+(1/31)];
% end
% for nContrast = 33:64
%     map(nContrast,:) = [map(nContrast-1,1)-(1/33), map(nContrast-1,2), map(nContrast-1,3)-(1/33)];
% end
% 
%             
% nBin=6;
% iSession=1;
% 
% 
% %----------------------------------------- Figure 1 ----------------------------------------------------%
% % scatter plots: probability of participation as a function of benefits and costs + indifference curves %
% Weight=figure;set(Weight, 'color', 'white');
% Weight.Units = 'centimeters'; Weight.Position = [ 1 1 40 20];
% 
% 
% 
%   % data names
%    taskNames = {'RE','PE'};
%    % taskNames = {'RE','PE','RP'};
%    benefitNames = {'récompenses','punitions','récompenses'};
%    costNames = {'efforts','efforts','punitions'};
%    tradeoffNames = {'récompenses/efforts','punitions/efforts','récompenses/punitions'};   
%    dimensionNames = {'récompenses','punitions','effort'};   
% 
%  
% for iTask = 1
%     
%        subplot(2,3,iTask)
%        hold on;box on;
%        % data
%        for iGroup = 1:2
%             eval([ 'benefit' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.weight.weight' taskNames{iTask} '.descriptiveStat.benefitStandardValue,1) ;' ]);
%             eval([ 'cost' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.weight.weight' taskNames{iTask} '.descriptiveStat.costStandardValue,1) ;' ]);
%             eval([ 'choice' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.weight.weight' taskNames{iTask} '.descriptiveStat.mean_goChoice_benefitStandardValue  ,1) ;' ]);
%             eval([ 'errorChoice' num2str(iGroup) ' = sem(groupResult.(groupNames{groupCompare(iGroup)}).group.weight.weight' taskNames{iTask} '.descriptiveStat.mean_goChoice_benefitStandardValue ,1) ;' ]);
%             eval([ 'prediction' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.weight.weight' taskNames{iTask} '.inferentialStat.mean_predicted_goChoice_benefitStandardValue,1) ;' ]);
% %             eval([ 'k' taskNames{iTask} ' = groupResult.(groupNames{groupCompare(iGroup)}).group.choice2D.statInferential.k' taskNames{iTask} '(1);' ]);
% %             eval([ 'temperature' taskNames{iTask} ' = groupResult.(groupNames{groupCompare(iGroup)}).group.choice2D.statInferential.temperature' taskNames{iTask} '(1);' ]);
% %             eval([ 'choiceFunction' iGroup ' = @(benefit,cost) 1./(1+exp(-( (k' taskNames{iTask} '*benefit-cost)./temperature' taskNames{iTask} '))); predictionH1 = [];']);
% 
%        end
%        
%        % plots
%        plot(benefit1 , prediction1,'b','Linewidth',2);
%        plot(benefit2 , prediction2,'r','Linewidth',2);
%        h1 = errorbar( benefit1 , choice1 , errorChoice1 ,'b','LineWidth',2,'LineStyle','none','Marker','o','MarkerFaceColor','auto');
%        h2 = errorbar( benefit2 , choice2 , errorChoice2 ,'r','LineWidth',2,'LineStyle','none','Marker','o','MarkerFaceColor','auto');
%  
%        % axis
%        axis = gca;
%        axis.TickLength = [0 0];
% %        axis.XTick = dimension; axis.XTickLabel = num2cell(round(dimension,1));
%        axis.YLim = [0 1]; % axis.XLim = [min(dimension)-0.1*(max(dimension)-min(dimension))  max(dimension)+0.1*(max(dimension)-min(dimension))];
%        xlabel(benefitNames{iTask},'FontSize',14);
%        ylabel('proportion d''acceptation(%)','FontSize',14);
%        title(tradeoffNames{iTask},'FontSize',14);
%  
%         % legend
%         legend([h1 h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
%        
%        % text
% end
% 
% 
% % store figure handles
% fig = Weight;
% 




end