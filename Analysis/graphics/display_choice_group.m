function [fig] = display_choice_2(dataTable,descriptiveTable,inferentialTable,option)
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



%% FIGURES %%
%__________________________________________________________________________
% Figure Specifications
title_size = 12;
axis_size=14;tick_size=10;
legend_size=14;
xSize= 40; ySize= 20;
x_text = 80; y_text = 80; y_shift = 10; text_color = '{0 0 0}';
star_size = 14;marker_size=0.3;
            
nBin=6;
iSession=1;


%----------------------------------------- Figure 1 ----------------------------------------------------%
% scatter plots & line plots: choice proportion as a function of delta value rating %
Choice=figure;set(Choice, 'color', 'white');
Choice.Units = 'centimeters'; Choice.Position = [ 1 1 40 20];

subTask = {'writtenReward','writtenPunishment','writtenEffort'};
subTaskAcronym = {'R','P','E'};
subTaskNames = {'récompenses','punitions','efforts'};


for iSubTask = [1,3]
    
    % data
    
    dataDV1 = nanmean( groupResult.(groupNames{groupCompare(1)}).group.choice.(subTask{iSubTask}).descriptiveStat.differenceZValue ,1)   ;
    dataChoice1 = nanmean( groupResult.(groupNames{groupCompare(1)}).group.choice.(subTask{iSubTask}).descriptiveStat.mean_sideChoice_Value  ,1 )  ;
    errorChoice1 = sem( groupResult.(groupNames{groupCompare(1)}).group.choice.(subTask{iSubTask}).descriptiveStat.mean_sideChoice_Value  ,1 )  ;
    dataAbsDV1 = nanmean( groupResult.(groupNames{groupCompare(1)}).group.choice.(subTask{iSubTask}).descriptiveStat.absoluteDifferenceZValue  ,1 )   ;
    dataRT1 = nanmean( groupResult.(groupNames{groupCompare(1)}).group.choice.(subTask{iSubTask}).descriptiveStat.mean_rt_absValue ,1)  ;
    errorRT1 = sem( groupResult.(groupNames{groupCompare(1)}).group.choice.(subTask{iSubTask}).descriptiveStat.mean_rt_absValue ,1)  ;
    
    dataDV2 = nanmean( groupResult.(groupNames{groupCompare(2)}).group.choice.(subTask{iSubTask}).descriptiveStat.differenceZValue ,1 )   ;
    dataChoice2 = nanmean( groupResult.(groupNames{groupCompare(2)}).group.choice.(subTask{iSubTask}).descriptiveStat.mean_sideChoice_Value ,1   )  ;
    errorChoice2 = sem( groupResult.(groupNames{groupCompare(2)}).group.choice.(subTask{iSubTask}).descriptiveStat.mean_sideChoice_Value ,1  )  ;
    dataAbsDV2 = nanmean( groupResult.(groupNames{groupCompare(2)}).group.choice.(subTask{iSubTask}).descriptiveStat.absoluteDifferenceZValue ,1  )   ;
    dataRT2 = nanmean( groupResult.(groupNames{groupCompare(2)}).group.choice.(subTask{iSubTask}).descriptiveStat.mean_rt_absValue ,1)  ;
    errorRT2 = sem( groupResult.(groupNames{groupCompare(2)}).group.choice.(subTask{iSubTask}).descriptiveStat.mean_rt_absValue ,1)  ;

 

    % subplot by Task
    subplot(2,3,iSubTask) % choices
    hold on;
    
    % parameters
    beta1 = nanmean(groupResult.(groupNames{groupCompare(1)}).group.choice.summary.(['temperature' subTaskAcronym{iSubTask} ])) ;
    beta2 = nanmean(groupResult.(groupNames{groupCompare(2)}).group.choice.summary.(['temperature' subTaskAcronym{iSubTask} ])) ;

%     BCA = texlabel(['BCA = ' num2str(round(result.(subTaskName{iSubTask}).inversion.fit.BCA ,2)) ]);
%     betaText =  texlabel(['beta_' subTaskAcronym{iSubTask} ' = ' num2str(round(beta,2))]);
    
    % model 
    f = @(x,beta) 1/(1+exp(-x/beta));
    x=[-100:0.1:100];
    for i=1:length(x)
        y1(i)=f(x(i),beta1);
        y2(i)=f(x(i),beta2);
    end
    plot(x,y1,'b','LineWidth',2);
    plot(x,y2,'r','LineWidth',2);
    
    % data
    h = scatter(dataDV1,dataChoice1,50,'b','filled');
    h2 = scatter(dataDV2,dataChoice2,50,'r','filled');
    errorbar(dataDV1,dataChoice1,errorChoice1,'b','LineStyle','none');
    errorbar(dataDV2,dataChoice2,errorChoice2,'r','LineStyle','none');

    xlabel('Valeur A - B (zscore)','FontSize',14);
    ylabel('Probabilité du choix A','FontSize',14);set(gca,'YTick',[0:0.1:1]); set(gca,'YTickLabel',{'0';'0.1';'0.2';'0.3';'0.4';'0.5';'0.6';'0.7';'0.8';'0.9';'1'});
%     text(-2.8,0.9,betaText,'FontSize',12,'Color','k','HorizontalAlignment','left');
%     text(-2.8,0.7,BCA,'FontSize',12,'HorizontalAlignment','left');
    plot([-10 10],[0 0],'k');
    plot([-10 10],[1 1],'k');
    % legend('Modèle','Choix','Données','Location',[0.65 0.3 0.1 0.1]);
    axis([-3 3 -0.2 1.2]);
    title(subTaskNames{iSubTask},'FontSize',14);
    
    % legend
    legend([h h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');

    
    subplot(2,3,iSubTask+3) % RT
    hold on;
    
    % parameters
    
    % model
    % data
    
    h = scatter(dataAbsDV1,dataRT1,20,'b','filled');
    h2 = scatter(dataAbsDV2,dataRT2,20,'r','filled');
    errorbar(dataAbsDV1,dataRT1,errorRT1,'b');
    errorbar(dataAbsDV2,dataRT2,errorRT2,'r');
    
    xlabel(' Valeur A - B (zscore)','FontSize',14);
    ylabel('temps de réaction (ms.)','FontSize',14);
    
    


end

% store figure handles
fig = Choice;





end