function [fig] = display_grip_2(groupNames,groupCompare,groupResult,groupData)
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



%----------------------------------------- Figure 1 ----------------------------------------------------%
% scatter plots: probability of participation as a function of benefits and costs + indifference curves %
Grip=figure;set(Grip, 'color', 'white');
Grip.Units = 'centimeters'; Grip.Position = [ 1 1 30 14];



    
       subplot(1,2,1)
       hold on;box on;
       % data
       for iGroup = 1:2
            eval([ 'data' num2str(iGroup) ' = groupResult.(groupNames{groupCompare(iGroup)}).group.grip.interDimensions.descriptiveStat.mean_normalizedForcePeak_incentive_Gain   ;' ]);
            eval([ 'data' num2str(iGroup) '(data' num2str(iGroup) '==0) = NaN   ;' ]);
            eval([ 'meanData' num2str(iGroup) ' = nanmean(data' num2str(iGroup) ',1)   ;' ]);
            eval([ 'errorData' num2str(iGroup) ' = sem(data' num2str(iGroup) ',1)   ;' ]);
            eval([ 'prediction' num2str(iGroup) ' = groupResult.(groupNames{groupCompare(iGroup)}).group.grip.interDimensions.descriptiveStat.mean_normalizedForcePeak_incentive_Gain   ;' ]);
            eval([ 'prediction' num2str(iGroup) '(data' num2str(iGroup) '==0) = NaN   ;' ]);
            eval([ 'meanPrediction' num2str(iGroup) ' = nanmean(data' num2str(iGroup) ',1)   ;' ]);       end
       
       % plots
       h1 = errorbar(  meanData1 , errorData1 ,'b','LineWidth',2,'LineStyle','-','Marker','o','MarkerFaceColor','auto');
       h2 = errorbar(  meanData2 , errorData2 ,'r','LineWidth',2,'LineStyle','-','Marker','o','MarkerFaceColor','auto');
 
       % axis
       axis = gca;
       axis.TickLength = [0 0];
       axis.XTick = [1:6]; axis.XTickLabel = {'0.01';'0.10';'0.50';'1';'5';'20'};
%        axis.YLim = [0 1]; % axis.XLim = [min(dimension)-0.1*(max(dimension)-min(dimension))  max(dimension)+0.1*(max(dimension)-min(dimension))];
       xlabel('Gain  (€)','FontSize',14);
       ylabel('Force (%. calibration)','FontSize',14);
 
        % legend
        legend([h1 h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
       
%               
% subplot(1,2,2)
%        hold on;box on;
%        % data
%        for iGroup = 1:2
%             eval([ 'data' num2str(iGroup) ' = groupResult.(groupNames{groupCompare(iGroup)}).group.grip.interDimensions.descriptiveStat.mean_normalizedForcePeak_incentive_Loss   ;' ]);
%             eval([ 'data' num2str(iGroup) '(data' num2str(iGroup) '==0) = NaN   ;' ]);
%             eval([ 'meanData' num2str(iGroup) ' = nanmean(data' num2str(iGroup) ',1)   ;' ]);
%             eval([ 'errorData' num2str(iGroup) ' = sem(data' num2str(iGroup) ',1)   ;' ]);
%        end
%        
%        % plots
%        h1 = errorbar(  meanData1 , errorData1 ,'b','LineWidth',2,'LineStyle','-','Marker','o','MarkerFaceColor','auto');
%        h2 = errorbar(  meanData2 , errorData2 ,'r','LineWidth',2,'LineStyle','-','Marker','o','MarkerFaceColor','auto');
%  
%        % axis
%        axis = gca;
%        axis.TickLength = [0 0];
%        axis.XTick = [1:6]; axis.XTickLabel = {'0.01';'0.10';'0.50';'1';'5';'20'};
% %        axis.YLim = [0 1]; % axis.XLim = [min(dimension)-0.1*(max(dimension)-min(dimension))  max(dimension)+0.1*(max(dimension)-min(dimension))];
%        xlabel('Pertes  (€)','FontSize',14);
%        ylabel('Force (%. calibration)','FontSize',14);
%  
%         % legend
%         legend([h1 h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
        
% store figure handles
fig = Grip;





end