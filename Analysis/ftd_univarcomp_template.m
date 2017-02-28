%% ftd_bivarcomp_template

%%% 3.1 statistics
% variable definition
ndim = 7; ngroup = 2;
nbin = [ ndim ,2, ngroup , nsub ];
taskname = 'rating';
Y = nan(nbin);

varnames = {'group','mean_r','mean_e',...
            'mean_food','mean_good','mean_motor','mean_cognitive',...
            'extrem_r','extrem_e',...
            'responsetime_r','responsetime_e'};
nvar = numel(varnames);
stat.(taskname) = array2table(nan(nsub,nvar)) ;
stat.(taskname).Properties.VariableNames = varnames;

% first-level statistics
for isub = 1:nsub
    try
    % select
    tab = data{isub}.battery.table; 
    selection =  (tab.task==taskname) ;
    tab = tab(selection,:);

    % variables
    rating = tab.rating/100 ;
    dim = tab.dimension;
    type = tab.itemSubtype;
    rt = tab.reactionTime;

    % stats
    [~,subdim] = ismember(unique(dim),dimensionList);
    igroup = (isub>n_control) + 1;
    stat.(taskname).group(isub) = igroup;
        % mean
        ysub = nan(1,ndim);
        ysub(subdim) = tools.tapply(rating,{dim},@nanmean);
        stat.(taskname){isub,[2 3]} = ysub([1 4]);
        ysub = nan(1,4);
        ysub = tools.tapply(rating(dim==dimensionList(1) | dim==dimensionList(4)),{type(dim==dimensionList(1) | dim==dimensionList(4))},@nanmean);
        stat.(taskname){isub,[4:7]} = ysub;
        % extremes
        ysub = nan(1,ndim);
        for idim = subdim'
            subrating = rating(dim==dimensionList(idim));
            % method 1
%             rge = range(subrating);
%             rank = quantileranks(subrating,3);
%             freq = histcounts(subrating,3,'Normalization','Probability');
%             extreme = freq(1) + freq(3);
            % method 2
            count = (subrating==min(subrating) | subrating==max(subrating));
            extreme = mean(count);
            
            ysub(idim) =  extreme;
        end
        stat.(taskname){isub,[8 9]} = ysub([1 4]);
        % responsetime
        ysub = nan(1,ndim);
        ysub(subdim) = tools.tapply(rt,{dim},@nanmean);
        stat.(taskname){isub,[10 11]} = ysub([1 4]);

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end




% check
disp(stat.rating);

%%% 3.2 display
% parameters definition
fig = figure; set(fig,'Name','rating_ftd');
hold on;clear h;
yvar = {'mean_food','mean_good','mean_motor','mean_cognitive'};
nvar = numel(yvar);
xpos = [1:4];
xticktext = {'food','good','cognitive','motor'};
ytext = 'mean rating (%)';

% plot
for ivar=1:nvar
    for ig = [1 2]
        xx = xpos(ivar) +(ig-1)*0.3;
        y = stat.(taskname).(yvar{ivar})(stat.(taskname).group==ig);
        yy = mean(y);
        zz = sem(y);
        % bar plot
%         [ h(ig)  ] = barplot( xx ,yy,zz, col{ig} );
%         h(ig).BarWidth = 0.3;
        % box plot
        boxplot(y,'boxstyle','outline','colors',col{ig},'medianstyle','line','symbol','','whisker',1,'widths',0.3,'positions',xx);
        if ivar==1; b = findobj('Tag','Box','-and','Color',col{ig}); h(ig)=b(1) ; end
        % dot plot
        dotdensity(xx*ones(numel(y),1)',y','dotEdgeColor',col{ig},'dotFaceColor','none','dotSize',5);
    end
    % inferential stats
    y1=stat.(taskname).(yvar{ivar})(stat.(taskname).group==1);
    y2=stat.(taskname).(yvar{ivar})(stat.(taskname).group==2);
    [h1,p] = ttest2(y1,y2);
    pos = [xpos(ivar) , xpos(ivar) + 0.3];
    [s] = sigstar(pos,p);
    
end

% legending
legend([h(1) h(2)],cellstr(groupList));
ax = gca; 
ax.TickLength = [0 0];
ax.XTick = xpos +0.15 ;
ax.XTickLabel = xticktext;
ax.XLim = [0 xpos(end)+1];       
ax.YLim = [0 1.25];       
ylabel(ytext); 
        
% format
setFigProper('FontSize',20,'LineWidth',2);