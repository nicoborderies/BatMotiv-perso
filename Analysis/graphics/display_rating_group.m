function [fig] = display_rating_2(dataTable,descriptiveTable,inferentialTable,option)
%% Specifications
    subTask = {'writtenReward','writtenPunishment','writtenEffort'};
    subTaskName = {'reward(text)','punishment(text)','effort(text)'};

    subTaskSubCategories = {'alimentaires','~alimentaires';
                            '~sensoriels','sensoriels';
                            'cognitifs','moteurs'};


%% Figure parameters
    colorTreatment =[0, 0 ,0 ; 1, 0, 1];
    colorSession = [0, 0 ,0 ; 0, 0, 1];
    title_size = 12;
    axis_size=14;tick_size=10;
    legend_size=14;
    xSize= 40; ySize= 20;
    x_text = 80; y_text = 80; y_shift = 10; text_color = '{0 0 0}';
    star_size = 14;marker_size=0.3;


%% Fig 1 : Histograms of ratings
%--------------------------------
% define fig handle
    ratingFig = figure;set(ratingFig, 'color', 'white');
    ratingFig.Units = 'centimeters'; ratingFig.Position = [ 1 1 30 15];

% iteration across subtask
    for iSubTask = 1:numel(subTask)
        
        % plot selection
            subplot(1,3,iSubTask);
            hold on;

        % data to plot
            design = groupData.(groupNames{groupCompare(1)}).group.rating.table;
            design = design(ismember(design.dimension,subTask{iSubTask}),:);
            data1 = design.rating  ;

            design = groupData.(groupNames{groupCompare(2)}).group.rating.table;
            design = design(ismember(design.dimension,subTask{iSubTask}),:);
            data2 = design.rating  ;

            h = histogram(data1,20,'Normalization','probability','FaceColor','b','FaceAlpha',0.5);
            h2 = histogram(data2,20,'Normalization','probability','FaceColor','r','FaceAlpha',0.5);

            % [frating,x] = ksdensity(rating);
            % plot(x,frating*max(histcounts(rating,10)/10)/max(frating),'r');

        % legend
            xlabel('ratings','FontSize',14);
            ylabel('proportion (%)','FontSize',14);
            title(subTask{iSubTask},'FontSize',14);

            legend([h h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');

        % axis 
            ax = get(gca);box on;
            set(gca,'YLim',[0 1],...
            'XLim',[-10 110],...
            'XTick',[0 50 100],...
            'XTickLabel',{'0';'50';'100'},...
            'TickLength',[0 0]);

    end

%% Save fig handle
    fig{1} = ratingFig;

end

