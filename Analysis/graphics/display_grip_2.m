function [fig] = display_grip_2(dataTable,descriptiveTable,inferentialTable,option)
%% Specifications
    valenceAcronym = {'R','P'};
    valenceName = {'Gain','Loss'};
    
    if nargin>3
        if isfield(option,'withinSubFactor')
            withinSubFactor = option.withinSubFactor;
            if isequal(option.withinSubFactor,{'constant'});
                dataTable.constant = ones(height(dataTable),1); 
            end
        else
            withinSubFactor =  {'sessionNumber'};
        end
    end
    

%% Figure parameters

color ={[0, 0 ,0] ; [1, 0, 1]};
colorSession = [0, 0 ,0 ; 0, 0, 1];
colorGroup = {'b','r'};
title_size = 12;
axis_size=14;tick_size=10;
legend_size=14;
xSize= 40; ySize= 20;
x_text = 80; y_text = 80; y_shift = 10; text_color = '{0 0 0}';
star_size = 14;marker_size=0.3;



%% Fig 1: force dependency to incentive
%--------------------------------

% define fig handle
gripFig1 =figure;set(gripFig1, 'color', 'white');
gripFig1.Units = 'centimeters'; gripFig1.Position = [ 1 1 20 15];


% iteration across withinSubFactor
    factor = dataTable.(withinSubFactor{1});
    levels = unique(factor); if ~iscell(levels);levels = num2cell(levels);end
        for iLevels = 1:numel(levels)

% iteration across valence
    i=0;
    for iValence = [1,-1]
        i = i+1;
        % plot selection
        % -------------- errorbar of choice|benefits
            subplot(1,2,i);  hold on;
            % data to plot
                KC = nanmean(inferentialTable.(['grip_kE' ]));
                stdKC = nanstd(inferentialTable.(['grip_kE' ]));
            
                selection =  ismember(dataTable.group,{'CONTROL','CITALOPRAM'})...
                                & factor == levels{iLevels}...
                                & dataTable.incentiveSign==iValence;  % select group
                subjects = dataTable.subject; subList = unique(subjects(selection));
                INCENTIVE = []; FORCE=[];PREDICTION = [];
                for iSub = 1:numel(subList)
                     incentive = dataTable.incentiveValue(selection & subjects==subList(iSub));   
                     force = dataTable.normalizedForcePeak(selection & subjects==subList(iSub));   
                     prediction = dataTable.predicted_forcePeak(selection & subjects==subList(iSub))...
                         ./dataTable.maxObservedForcePeak(selection & subjects==subList(iSub)) ;  
                     FORCE(iSub,:) = tools.tapply(force,{incentive},@nanmean);
                     PREDICTION(iSub,:) = tools.tapply(prediction,{incentive},@nanmean);
                end
                K = nanmean(inferentialTable.(['grip_k' valenceAcronym{i}]));
                stdK = nanstd(inferentialTable.(['grip_k' valenceAcronym{i} ]));

                
            % plot
                h = errorbar( 1:numel(unique(incentive)) , nanmean(FORCE,1) , sem(FORCE,1) );                
                scatter( 1:numel(unique(incentive)) , nanmean(FORCE,1) , 'filled','MarkerFaceColor',color{iLevels});
                plot( 1:numel(unique(incentive)) , nanmean(PREDICTION,1) , '--','Color',color{iLevels} );
                h.LineStyle = 'none';  h.LineWidth = 1.5; h.Color = [0 0 0]; 
                
            % legend
                xx = xlim; xx = 1; yy = 0.9;
                xlabel(texlabel([ ' incentive (€) ' ]),'FontSize',14);
                ylabel('force (% MVC)','FontSize',14);
                title(valenceName{i},'FontSize',14);
                text(xx,yy,texlabel(['k_' valenceAcronym{i} ' = ' num2str(round(K,2)) ' (+/- ' num2str(round(stdK,2)) ')']));
                text(xx,yy-0.1,texlabel(['k_E = ' num2str(round(KC,2)) ' (+/- ' num2str(round(stdKC,2)) ')']));

                %             legend([h h2],{ groupNames{groupCompare(1)} , groupNames{groupCompare(2)} },'Location','northoutside');
            % axis 
                ax = get(gca);box off; ax.FontName = 'FixedWidth';ax.LabelFontSizeMultiplier = 0.9;
                xticks = unique(incentive);
                set(gca,...
                    'XTick',[1:6],...
                    'XTickLabel',cellfun(@num2str,num2cell(xticks),'UniformOutput',0),...
                    'YLim',[0 1],...
                    'TickLength',[0 0]);
          end
    end

    
%% Save fig handle
    fig{1} = gripFig1; 





end