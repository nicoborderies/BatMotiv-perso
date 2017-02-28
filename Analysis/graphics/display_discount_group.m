function [fig] = display_discount_group(groupNames,groupCompare,groupResult,groupData)
%% INITIALIZE 
    % variable definition
    if nargin <3;
    fig = [];end

    subTaskName={'writtenReward', 'writtenPunishment', 'writtenEffort'}; % to use as field name in data
    submanipAcronym = {'R','P','E'};
    dimNames = {'delta(reward)','delta(punishment)','delta(effort)'};
    
    
    
    % task selection
    nSubTask = length(subTaskName);
    selectSubTask=[];
    selectSubTask=[1];
%     for iSubTask = 1:nSubTask
%         if isfield(data,subTaskName{iSubTask})
%             if ~isempty(data.(subTaskName{iSubTask}))
%                 selectSubTask=[selectSubTask iSubTask];
%             end
%         end
%     end

    % session selection
%     if nargin<3;
        sessionList = []; 
%     end
    

%% FIGURES %%
%__________________________________________________________________________
% Figure Specifications
colorTreatment =[0, 0 ,0 ; 1, 0, 1];
colorSession = [0, 0 ,0 ; 0, 0, 1];
title_size = 12;
axis_size=14;tick_size=10;
legend_size=14;
xSize= 40; ySize= 20;
x_text = 80; y_text = 80; y_shift = 10; text_color = '{0 0 0}';
star_size = 14;marker_size=0.3;


map(1,:) = [1 0 0]; % colormap
for nContrast = 2:32
    map(nContrast,:) = [map(nContrast-1,1), map(nContrast-1,2)+(1/31), map(nContrast-1,3)+(1/31)];
end
for nContrast = 33:64
    map(nContrast,:) = [map(nContrast-1,1)-(1/33), map(nContrast-1,2), map(nContrast-1,3)-(1/33)];
end


%----------------------------------------- Figure 1 ----------------------------------------------------%
% scatter plots: probability of selecting delayed option as a function of reward and delay + indifference curves %
Discount=figure;set(Discount, 'color', 'white');   
Discount.Units = 'centimeters'; Discount.Position = [ 1 1 30 14];


for iTask = selectSubTask  % iteration across tasks
      

   %text
%     Kt = result.(Task{nTask}).inversion.parameters{1, 1}.Kt  ;
%     temperature = result.(Task{nTask}).inversion.parameters{1, 1}.temperature  ;
%     switch result.(Task{nTask}).inversion.specifications.discount
%         case 'quasiHyperbolic'
%             gamma = result.(Task{nTask}).inversion.parameters{1, 1}.gamma  ;
%     end
% BCA = texlabel(['BCA = ' num2str(round(result.(Task{nTask}).inversion.fit.BCA,2)) ]);
                
      
% ----------------------------------------- Figure 2 ----------------------------------------------------%
% errorbar plots: probability of selecting delayed option as a function of reward / delay  %

        
       subplot(2,3,iTask)
       hold on;box on;
       % data
      for iGroup = 1:2
           eval([ 'dimension' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.discount.' subTaskName{iTask} '.descriptiveStat.differenceZValue,1) ;' ]);
           eval([ 'perf' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.discount.' subTaskName{iTask} '.descriptiveStat.mean_delayedChoice_differenceZValue,1)   ;' ]);
           eval([ 'errorPerf' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.discount.' subTaskName{iTask} '.descriptiveStat.error_delayedChoice_differenceZValue,1)   ;' ]);
           eval([ 'prediction' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.discount.' subTaskName{iTask} '.inferentialStat.mean_predicted_delayedItemChoice_differenceZValue,1) ;' ]);
      end
      
       % plots
       h1 = plot(dimension1 , prediction1,'b','LineWidth',2);
       hError = errorbar( dimension1 , perf1 , errorPerf1 ,'b','LineWidth',2);
       hError.LineStyle = 'none';
       hError.Marker ='o';hError.MarkerFaceColor = 'auto';
       h2 = plot(dimension2 , prediction2,'r','LineWidth',2);
       hError = errorbar( dimension2 , perf2 , errorPerf2 ,'r','LineWidth',2);
       hError.LineStyle = 'none';
       hError.Marker ='o';hError.MarkerFaceColor = 'auto';
       
       % axis
       axis = gca;
       axis.TickLength = [0 0];
       axis.YLim = [0 1]; %axis.XLim = [0 100];
%        axis.XLim = [min(dimension1)  max(dimension)];
       xlabel(dimNames{1,iTask},'FontSize',14);
       ylabel('choix de l''option différée(%)','FontSize',12);
       
        % legend
        legend([h1 h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
       
       
       % text
%        text_temperature =  texlabel(['beta_' parametersDim{nTask} ' = ' num2str(round(temperature,2))]);
%        text(0*(max(dimension)-min(dimension))+min(dimension),0.9,text_temperature,'FontSize',14);
%        BCA = texlabel(['BCA = ' num2str(round(result.(Task{nTask}).inversion.fit.BCA,2)) ]);
%        text(0*(max(dimension)-min(dimension))+min(dimension),0.75,BCA,'FontSize',14);
               

       subplot(2,3,iTask+3)
       hold on;box on;
       % data
         for iGroup = 1:2
           eval([ 'dimension' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.discount.' subTaskName{iTask} '.descriptiveStat.delay,1) ;' ]);
           eval([ 'perf' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.discount.' subTaskName{iTask} '.descriptiveStat.mean_delayedChoice_delay,1)   ;' ]);
           eval([ 'errorPerf' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.discount.' subTaskName{iTask} '.descriptiveStat.error_delayedChoice_delay,1)   ;' ]);
           eval([ 'prediction' num2str(iGroup) ' = nanmean(groupResult.(groupNames{groupCompare(iGroup)}).group.discount.' subTaskName{iTask} '.inferentialStat.mean_predicted_delayedItemChoice_delay,1) ;' ]);
         end
         
       % plots
       h1 = plot(dimension1 , prediction1,'b','LineWidth',2);
       hError = errorbar( dimension1 , perf1 , errorPerf1 ,'b','LineWidth',2);
       hError.LineStyle = 'none';
       hError.Marker ='o';hError.MarkerFaceColor = 'auto';
       h2 = plot(dimension2 , prediction2,'r','LineWidth',2);
       hError = errorbar( dimension2 , perf2 , errorPerf2 ,'r','LineWidth',2);
       hError.LineStyle = 'none';
       hError.Marker ='o';hError.MarkerFaceColor = 'auto';
       
       % axis
       axis = gca;
       axis.XScale = 'log';
       axis.TickLength = [0 0];
%        dimRange = round(unique(dimension));
%        axis.XTick = dimRange(1:2:end); 
%        axis.XTickLabel = num2cell(dimRange(1:2:end));
       axis.YLim = [0 1];axis.XLim = [0 1000];
       xlabel('delay','FontSize',14);
       ylabel('choix de l''option différée(%)','FontSize',12);
       
        % legend
        legend([h1 h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
       
       
       %text
%        switch (1000/Kt)>365
%            case 0
%                text_Kt =  texlabel(['tau_' parametersDim{nTask} ' = ' num2str(round(1000/Kt,1)) ' jours']);
%            case 1
%                text_Kt =  texlabel(['tau_' parametersDim{nTask} ' = ' num2str(round(1000/(Kt*365.25),1)) ' ans']);
%        end       
%        text_gamma =  texlabel(['gamma_' parametersDim{nTask} ' = ' num2str(round((1/gamma),2)) ' ' parametersDim{nTask} '.']);
%        text(10,0.9,text_Kt,'FontSize',14);
%        text(10,0.8,text_gamma,'FontSize',14);
    
end

% store figure handles
fig{1} = Discount;
 
end

