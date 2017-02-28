%% ftd_bivarcomp_template

%%% .1 statistics
% variable definition
ndim = 4; ngroup = 2;
taskname = 'choice';
nqtle = 6;
nbin = [ nqtle,ndim,nsub ];
GROUP = nan(nsub,1);
X = nan(nbin);
Y = nan(nbin);
Y2 = nan(nbin);

% first-level statistics
for isub = 1:nsub
    try
    % select
    tab = data{isub}.battery.table; 
    selection =  (tab.task==taskname) ;
    tab = tab(selection,:);

    % variables
    dv = tab.differenceRatingValue/100 ;
    dim = tab.dimension;
    choice = (tab.sideChoice+1)/2;
    choice2 = tab.predicted_sideChoice;

    % stats
    [~,subdim] = ismember(unique(dim),dimensionList);
    igroup = (isub>n_control) + 1;
    GROUP(isub)=igroup;
    % histogram
    ysub = nan(1,ndim);
    for idim = subdim'
        select = dim==dimensionList(idim);
        xbin = quantileranks(dv(select),nqtle);
        X(:,idim,isub) = tools.tapply(dv(select),{xbin},@nanmean);
        Y(:,idim,isub)  = tools.tapply(choice(select),{xbin},@nanmean);
        Y2(:,idim,isub)  = tools.tapply(choice2(select),{xbin},@nanmean);
    end

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end

% check
disp(Y);

%%% .2 display
% parameters definition
fig = figure; set(fig,'Name','rating_ftd');
hold on;clear h;

% xticktext = {};
xtext = 'delta value (%)';
ytext = 'choice = right (%)';
xlimits = [0 1];
ylimits = [0 1];
dsub=3;

i=0;
for idim = [1 4]
    i=i+1;
    subplot(1,2,i);hold on;

    % plot
    for ivar=1:nvar
        for ig = [1 2]
            % data2plot
            xx = mean(X(:,idim,GROUP==ig),dsub);
            yy = mean(Y(:,idim,GROUP==ig),dsub);
            yy2 = mean(Y2(:,idim,GROUP==ig),dsub);
            zz = sem(Y2(:,idim,GROUP==ig),dsub);
            xx = xx +(ig-1)*0;

            % errobar plot
            [~,h(ig),l] = errorscat( xx ,yy, zz,col{ig});
            l.LineStyle='none';
            
            h(ig) = plot( xx ,yy2,'Color',col{ig});
        end
    end

    % legending
    legend([h(1) h(2)],cellstr(groupList));
    ax = gca; 
    % ax.TickLength = [0 0];
    % ax.XTick = [0:0.1:1] ;
    ax.XLim = xlimits;       
    ax.YLim = ylimits;     
    xlabel(xtext); 
    ylabel(ytext); 
    title(titletext{i}); 
end
% format
setFigProper('FontSize',20,'LineWidth',2);
        
% format
setFigProper('FontSize',20,'LineWidth',2);