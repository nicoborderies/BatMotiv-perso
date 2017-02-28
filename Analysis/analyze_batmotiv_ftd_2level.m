%% analyze_batmotiv_ftd_2level
%
% this script execute a group analysis of FTD/CONTROLS dataset with
% statistical comparions & visualization
%
% see subscripts:
%   - ftd_univarcomp_template
%
% Nicolas Borderies
% 01/2017

%% 1/ import
dataset = 'batmotiv_analysis_FTD_all_31_1_2017.mat';
load(dataset);

%% 2/ define parameters

% structure dataset 
data = [ groupData.CONTROL.subject , groupData.bvFTD.subject ];
result = [ groupResult.CONTROL.subject , groupResult.bvFTD.subject ];
nsub = numel(data);
n_control = 20;
n_ftd = 19;

% label variables
dimensionList = nominal({'writtenReward','visualReward','writtenPunishment','writtenEffort',...
                            'monetaryReward','monetaryPunishment','gripEffort',...
                            'RewardEffort','PunishmentEffort','RewardPunishment',...
                            'Gain','Loss','GainLoss','GainEmotion'});
taskList = nominal( {'rating','choice','weight','discount',...
                    'grip','gripIAPS','gripAccu','mental','learning'});
groupList = nominal({'CONTROL','FTD'});
group = [ repmat(groupList(1),1,n_control) , repmat(groupList(2),1,n_ftd) ];


% graphical parameters
col = { [1 1 1]*0.5 , [1 0 0]*1 };
    
%% 3/ rating analysis
%%% 3.1 averages
%%% 3.1.1 statistics
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

%%% 3.1.2 display
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

%% 3.2 histogram
%%% .1 statistics
% variable definition
ndim = 7; ngroup = 2;
taskname = 'rating';
nqtle = 10;
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
    rating = tab.rating/100 ;
    dim = tab.dimension;
    type = tab.itemSubtype;
    rt = tab.reactionTime;

    % stats
    [~,subdim] = ismember(unique(dim),dimensionList);
    igroup = (isub>n_control) + 1;
    GROUP(isub)=igroup;
    % histogram
    ysub = nan(1,ndim);
    for idim = subdim'
        subrating = rating(dim==dimensionList(idim));
        subprediction = tab.predicted_rating(dim==dimensionList(idim));
        if idim==4; subprediction = -subprediction;end
        interval=[0:0.1:1];
        freq = histcounts(subrating,[0:0.1:1],'Normalization','Probability');
        X(:,idim,isub) = interval(1:end-1)+0.05;
        Y(:,idim,isub) = freq;
        freq = histcounts(subprediction,[0:0.1:1],'Normalization','Probability');
        Y2(:,idim,isub) = freq;
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
fig = figure; set(fig,'Name','rating_ftd_2');
hold on;clear h;

% xticktext = {};
xtext = 'rating (%)';
ytext = 'frequency (%)';
titletext = {'reward','effort'};
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
            zz = sem(Y(:,idim,GROUP==ig),dsub);
            xx = xx +(ig-1)*0.03;
%             xx = xx +(ig-1)*0;

            % errobar plot
    %         [~,~,h(ig)] = errorscat( xx ,yy2, zz,col{ig});
            h(ig) = plot( xx ,yy2,'Color',col{ig});

            % histogram plot
    %         [ h(ig)  ] = bar( xx ,yy, 'EdgeColor',col{ig},'FaceColor',col{ig} );
            [ h(ig)  ] = bar( xx ,yy,0.3,'EdgeColor',col{ig},'FaceColor','none');

        end
    end

    % legending
    legend([h(1) h(2)],cellstr(groupList));
    ax = gca; 
    ax.TickLength = [0 0];
    ax.XTick = [0:0.1:1] ;
    ax.XLim = [0 1];       
    ax.YLim = [0 0.5];
    xlabel(xtext); 
    ylabel(ytext); 
    title(titletext{i}); 
end
        
% format
setFigProper('FontSize',20,'LineWidth',2);

%% 4/ choice-1D analysis
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
xtext = '\Delta value (%)';
ytext = 'choice = right (%)';
titletext = {'reward','effort'};
% xlimits = [0 1];
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
            zz = sem(Y(:,idim,GROUP==ig),dsub);
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
    % ax.XLim = [0 1];       
    % ax.YLim = [0 1.25];       
    xlabel(xtext); 
    ylabel(ytext); 
    title(titletext{i}); 
end
% format
setFigProper('FontSize',20,'LineWidth',2);
    
%% 5/ choice-2D analysis

%%%% benefit dependency
%%% .1 statistics
% variable definition
ndim = 4; ngroup = 2;
taskname = 'weight';
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
    value = tab.benefitItemRatingValue/100 ;
    accept = tab.isGoChoice;
    accept2 = tab.predicted_goChoice;


    % stats
    [~,subdim] = ismember(unique(dim),dimensionList);
    igroup = (isub>n_control) + 1;
    GROUP(isub)=igroup;
    select = ~isnan(value);
    xbin = quantileranks(value(select),nqtle);
    X(:,idim,isub) = tools.tapply(value(select),{xbin},@nanmean);
    Y(:,idim,isub)  = tools.tapply(accept(select),{xbin},@nanmean);
    Y2(:,idim,isub)  = tools.tapply(accept2(select),{xbin},@nanmean);

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end

% check
disp(Y);

%%% .2 display
% parameters definition
fig = figure; set(fig,'Name','weight_ftd');
hold on;clear h;

% xticktext = {};
xtext = 'reward value (%)';
ytext = 'choice = accept (%)';
% xlimits = [0 1];
ylimits = [0 1];
dsub=3;

i=1;
subplot(1,2,i);hold on;

% plot
for ivar=1:nvar
    for ig = [1 2]
        % data2plot
        xx = mean(X(:,idim,GROUP==ig),dsub);
        yy = mean(Y(:,idim,GROUP==ig),dsub);
        yy2 = mean(Y2(:,idim,GROUP==ig),dsub);
        zz = sem(Y(:,idim,GROUP==ig),dsub);
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
% ax.XLim = [0 1];       
ax.YLim = [0 1];    
xlabel(xtext); 
ylabel(ytext); 

%%%% cost dependency
%%% .1 statistics
% variable definition
ndim = 4; ngroup = 2;
taskname = 'weight';
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
    value = tab.costItemRatingValue/100 ;
    accept = tab.isGoChoice;
    accept2 = tab.predicted_goChoice;


    % stats
    [~,subdim] = ismember(unique(dim),dimensionList);
    igroup = (isub>n_control) + 1;
    GROUP(isub)=igroup;
    select = ~isnan(value);
    xbin = quantileranks(value(select),nqtle);
    X(:,idim,isub) = tools.tapply(value(select),{xbin},@nanmean);
    Y(:,idim,isub)  = tools.tapply(accept(select),{xbin},@nanmean);
    Y2(:,idim,isub)  = tools.tapply(accept2(select),{xbin},@nanmean);

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end

% check
disp(Y);

%%% .2 display
% parameters definition


% xticktext = {};
xtext = 'effort value (%)';
ytext = 'choice = accept (%)';
% xlimits = [0 1];
ylimits = [0 1];
dsub=3;

i=2;
subplot(1,2,i);hold on;

% plot
for ivar=1:nvar
    for ig = [1 2]
        % data2plot
        xx = mean(X(:,idim,GROUP==ig),dsub);
        yy = mean(Y(:,idim,GROUP==ig),dsub);
        yy2 = mean(Y2(:,idim,GROUP==ig),dsub);
        zz = sem(Y(:,idim,GROUP==ig),dsub);
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
% ax.XLim = [0 1];       
ax.YLim = [0 1];    
xlabel(xtext); 
ylabel(ytext); 

% format
setFigProper('FontSize',20,'LineWidth',2);

%% 6/ choice-IT analysis

%%%% benefit dependency
%%% .1 statistics
% variable definition
ndim = 4; ngroup = 2;
taskname = 'discount';
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
    v1 = tab.immediateItemRatingValue/100 ;
    v2 = tab.delayedItemRatingValue/100 ;
    dv = v2-v1;
    patient = tab.isDelayedItemChoice;
    patient2 = tab.predicted_delayedItemChoice;


    % stats
    [~,subdim] = ismember(unique(dim),dimensionList);
    igroup = (isub>n_control) + 1;
    GROUP(isub)=igroup;
    select = ~isnan(dv);
    xbin = quantileranks(dv(select),nqtle);
    X(:,idim,isub) = tools.tapply(dv(select),{xbin},@nanmean);
    Y(:,idim,isub)  = tools.tapply(patient(select),{xbin},@nanmean);
    Y2(:,idim,isub)  = tools.tapply(patient2(select),{xbin},@nanmean);

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end

% check
disp(Y);

%%% .2 display
% parameters definition
fig = figure; set(fig,'Name','discount_ftd');
hold on;clear h;

% xticktext = {};
xtext = ' \Delta reward value (%)';
ytext = 'choice = patient (%)';
% xlimits = [0 1];
ylimits = [0 1];
dsub=3;

i=1;
subplot(1,2,i);hold on;

% plot
for ivar=1:nvar
    for ig = [1 2]
        % data2plot
        xx = mean(X(:,idim,GROUP==ig),dsub);
        yy = mean(Y(:,idim,GROUP==ig),dsub);
        yy2 = mean(Y2(:,idim,GROUP==ig),dsub);
        zz = sem(Y(:,idim,GROUP==ig),dsub);
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
% ax.XLim = [0 1];       
ax.YLim = [0 1];    
xlabel(xtext); 
ylabel(ytext); 

%%%% cost dependency
%%% .1 statistics
% variable definition
ndim = 4; ngroup = 2;
taskname = 'discount';
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
    delay = tab.delay;
    patient = tab.isDelayedItemChoice;
    patient2 = tab.predicted_delayedItemChoice;


    % stats
    [~,subdim] = ismember(unique(dim),dimensionList);
    igroup = (isub>n_control) + 1;
    GROUP(isub)=igroup;
    select = ~isnan(delay);
    xbin = quantileranks(delay(select),nqtle);
    X(:,idim,isub) = tools.tapply(delay(select),{xbin},@nanmean);
    Y(:,idim,isub)  = tools.tapply(patient(select),{xbin},@nanmean);
    Y2(:,idim,isub)  = tools.tapply(patient2(select),{xbin},@nanmean);

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end

% check
disp(Y);

%%% .2 display
% parameters definition


% xticktext = {};
xtext = 'delay (days)';
ytext = 'choice = patient (%)';
% xlimits = [0 1];
ylimits = [0 1];
dsub=3;

i=2;
subplot(1,2,i);hold on;

% plot
for ivar=1:nvar
    for ig = [1 2]
        % data2plot
        xx = mean(X(:,idim,GROUP==ig),dsub);
        yy = mean(Y(:,idim,GROUP==ig),dsub);
        yy2 = mean(Y2(:,idim,GROUP==ig),dsub);
        zz = sem(Y(:,idim,GROUP==ig),dsub);
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
ax.XTick = [0 1 7 14 21 45 100 365 1000] ;
% ax.XLim = [0 1];       
ax.XScale = 'log';
ax.XLim = [0 1000];       
ax.YLim = [0 1];    
xlabel(xtext); 
ylabel(ytext); 

% format
setFigProper('FontSize',20,'LineWidth',2);


%% 7/ model-based preference analysis
%%% 3.1 statistics
% variable definition
ndim = 2; ngroup = 2;
nbin = [ ndim ,2, ngroup , nsub ];
taskname = 'battery';
Y = nan(nbin);

varnames = {'group','mu_R','mu_E','sd_R','sd_E','alpha','kD','bD','bR','bm','R2_rating','BCA_choice'};
nvar = numel(varnames);
stat.(taskname) = array2table(nan(nsub,nvar)) ;
stat.(taskname).Properties.VariableNames = varnames;

% first-level statistics
for isub = 1:nsub
    try
    % select
    tab = result{isub}.battery.inferential;

    % stats
    igroup = (isub>n_control) + 1;
    stat.(taskname).group(isub) = igroup;
    
    stat.(taskname){isub,[2:12]} = [ tab.muR tab.muE tab.stdR tab.stdE tab.alpha tab.kRD tab.bRD tab.bR tab.bm tab.R2(1) tab.BCA(2)];
    

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end




% check
disp(stat.battery);
groupstat = varfun(@nanmean,stat.battery,'GroupingVariables','group');
disp(groupstat);


%%% 3.2 display
% parameters definition
fig = figure; set(fig,'Name','grip_ftd');
hold on;clear h;
yvar = {'mu_R','mu_E','sd_R','sd_E','alpha','kD','bD','bR','bm','R2_rating','BCA_choice'};
nvar = numel(yvar);
xpos = [1:nvar];
xticktext = {'mu_R','mu_E','sd_R','sd_E','alpha','kD','bD','bR','bm','R2_rating','BCA_choice'};
ytext = 'parameters (au.)';

% plot
ip=[1:9];
subplot(1,4,[1 2 3]);hold on;
for ivar=ip
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
ax.XLim = [ min(xpos(ip))-0.5 max(xpos(ip))+0.5];       
ylabel(ytext);

% plot
ip=[10 11];
subplot(1,4,4);hold on;
for ivar=ip
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
ax.XLim = [ min(xpos(ip))-0.5 max(xpos(ip))+0.5];       
ylabel(ytext);


% format
setFigProper('FontSize',20,'LineWidth',2);

%% 8/ force analysis
%%%% 8.1 model-free behavioral results
%%% .1 statistics
% variable definition
ngroup = 2;
taskname = 'grip';
nqtle = 6;nqtle2 = 10;
nbin = [ nqtle,nsub ];
nbin2 = [ nqtle2,nsub ];
GROUP = nan(nsub,1);
X = nan(nbin);
Y = nan(nbin);
Yhat = nan(nbin);
X2 = nan(nbin2);
Y2 = nan(nbin2);
Yhat2 = nan(nbin2);
X3 = nan(nbin2);
Y3 = nan(nbin2);
Yhat3 = nan(nbin2);

% first-level statistics
for isub = 1:nsub
    try
    % select
    tab = data{isub}.battery.table; 
    selection =  (tab.task==taskname) ;
    tab = tab(selection,:);

    % variables
    value = tab.incentiveLevel;
    nt = tab.trialNumber;
    force = tab.normalizedForcePeak;
    force2 = tab.predicted_forcePeak./tab.maxObservedForcePeak;
    yank = tab.yankPeak./tab.maxObservedForcePeak;
    yank2 = tab.predicted_yankPeak./tab.maxObservedForcePeak;

    % stats
    igroup = (isub>n_control) + 1;
    GROUP(isub)=igroup;
    
    xbin = quantileranks(value,nqtle);
    X(:,isub) = tools.tapply(value,{xbin},@nanmean);
    Y(:,isub)  = tools.tapply(force,{xbin},@nanmean);
    Yhat(:,isub)  = tools.tapply(force2,{xbin},@nanmean);
    
    xbin = quantileranks(nt,nqtle2);
    X2(:,isub) = tools.tapply(nt,{xbin},@nanmean);
    Y2(:,isub)  = tools.tapply(force,{xbin},@nanmean);
    Yhat2(:,isub)  = tools.tapply(force2,{xbin},@nanmean);
    
    xbin = quantileranks(force,nqtle2);    xbin(xbin==0)=NaN;
    xbin2 = quantileranks(force2,nqtle2);    xbin2(xbin2==0)=NaN;
    X3(:,isub) = tools.tapply(force,{xbin},@nanmean);
    Y3(:,isub)  = tools.tapply(yank,{xbin},@nanmean);
    Yhat3(:,isub)  = tools.tapply(yank2,{xbin2},@nanmean);

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end

% check
% disp(Y);
% disp(Y2);
% disp(Y3);


%%% .2 display
% parameters definition
fig = figure; set(fig,'Name','grip_ftd');
hold on;clear h;


%%%% force-incentive plot
% xticktext = {};
xtext = ' incentive reward (€)';
ytext = 'force peak (%fmax)';
xlimits = [0 nqtle+1];
ylimits = [0 1];
dsub=2;
subplot(2,2,1);hold on;
% plot
for ivar=1:nvar
    for ig = [1 2]
        % data2plot
        xx = mean(X(:,GROUP==ig),dsub);
        yy = mean(Y(:,GROUP==ig),dsub);
        yy2 = mean(Yhat(:,GROUP==ig),dsub);
        zz = sem(Y(:,GROUP==ig),dsub);
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
ax.XTick = [1:nqtle] ;
ax.XTickLabel = {'0.01','0.2','0.50','1','5','20'} ;
ax.XLim = xlimits;       
ax.YLim = ylimits;    
xlabel(xtext); 
ylabel(ytext); 

%%%% force-ntrial plot
% xticktext = {};
xtext = ' trial number (n)';
ytext = 'force peak (%fmax)';
xlimits = [0 70];
ylimits = [0 1];
dsub=2;
subplot(2,2,2);hold on;
% plot
for ivar=1:nvar
    for ig = [1 2]
        % data2plot
        xx = mean(X2(:,GROUP==ig),dsub);
        yy = mean(Y2(:,GROUP==ig),dsub);
        yy2 = mean(Yhat2(:,GROUP==ig),dsub);
        zz = sem(Y2(:,GROUP==ig),dsub);
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
ax.XTick = [0:6:60] ;
ax.XLim = xlimits;       
ax.YLim = ylimits;    
xlabel(xtext); 
ylabel(ytext); 


%%%% yank-force plot
% xticktext = {};
xtext = 'force peak (%fmax)';
ytext = 'velocity peak (%fmax/s)';
xlimits = [0 1];
% ylimits = [0 1];
dsub=2;
subplot(2,2,4);hold on;
% plot
for ivar=1:nvar
    for ig = [1 2]
        % data2plot
        xx = mean(X3(:,GROUP==ig),dsub);
        yy = mean(Y3(:,GROUP==ig),dsub);
        yy2 = mean(Yhat3(:,GROUP==ig),dsub);
        zz = sem(Y3(:,GROUP==ig),dsub);
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
ax.XTick = [0:0.1:1] ;
ax.XLim = xlimits;       
% ax.YLim = ylimits;    
xlabel(xtext); 
ylabel(ytext); 


%%%% 8.2 model-free heart-rate response
%%% .1 statistics
% variable definition
ngroup = 2;
taskname = 'grip';
nqtle = 6;nqtle2 = 10;
nbin = [ nqtle,nsub ];
nbin2 = [ nqtle2,nsub ];
GROUP = nan(nsub,1);
X = nan(nbin);
Y = nan(nbin);
Yhat = nan(nbin);

% first-level statistics
for isub = 1:nsub
    try
    % select
    tab = result{isub}.grip.oxmodel.oxmatrix;
    estimate = result{isub}.grip.inferential;

    % variables
    mu = estimate.mu_bpm;
    value = tab.incentive;
    force = tab.force;
    hr = (tab.hr_effort - mu) ;

    % stats
    igroup = (isub>n_control) + 1;
    GROUP(isub)=igroup;
    
    % iterative regression method
    predictor = nanzscore([ones(size(hr)),value , force  ]);
    predictor(:,1) = 1;
    response = hr;
    % 1. regress response to force
    [beta,~,~,out] = nanglm(predictor(:,[1 3]),response,0);
    prediction = out.suffStat.gx;
    residuals = (response - prediction);
    response = residuals;
    % 2. regress response to incentive
    [beta2,~,~,out] = nanglm(predictor(:,[1 2]),response,0);
    hr2 = prediction + (out.suffStat.gx);
    
    xbin = quantileranks(value,nqtle);
    X(:,isub) = tools.tapply(value,{xbin},@nanmean);
    Y(:,isub)  = tools.tapply(hr,{xbin},@nanmean);
    Yhat(:,isub)  = tools.tapply(hr2,{xbin},@nanmean);
    

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end

% check
% disp(Y);
% disp(Y2);
% disp(Y3);


%%% .2 display
% parameters definition

%%%% force-incentive plot
% xticktext = {};
xtext = ' incentive reward (€)';
ytext = 'heart rate response (bpm)';
xlimits = [0 nqtle+1];
ylimits = [0 1];
dsub=2;
i=4;
subplot(2,2,3);hold on;
% plot
for ivar=1:nvar
    for ig = [1 2]
        % data2plot
        xx = mean(X(:,GROUP==ig),dsub);
        yy = mean(Y(:,GROUP==ig),dsub);
        yy2 = mean(Yhat(:,GROUP==ig),dsub);
        zz = sem(Y(:,GROUP==ig),dsub);
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
ax.XTick = [1:nqtle] ;
ax.XTickLabel = {'0.01','0.2','0.50','1','5','20'} ;
ax.XLim = xlimits;       
% ax.YLim = ylimits;    
xlabel(xtext); 
ylabel(ytext); 

% format
setFigProper('FontSize',20,'LineWidth',2);

%% 8.2 model-based analysis
%%% 3.1 statistics
% variable definition
ndim = 2; ngroup = 2;
nbin = [ ndim ,2, ngroup , nsub ];
taskname = 'grip';
Y = nan(nbin);

varnames = {'group','kr','ke','kf','tau','R2_force','R2_velocity'};
nvar = numel(varnames);
stat.(taskname) = array2table(nan(nsub,nvar)) ;
stat.(taskname).Properties.VariableNames = varnames;

% first-level statistics
for isub = 1:nsub
    try
    % select
    tab = result{isub}.grip.inferential;

    % stats
    igroup = (isub>n_control) + 1;
    stat.(taskname).group(isub) = igroup;
    
    stat.(taskname){isub,[2:7]} = [ tab.kR tab.kE tab.kF tab.tau tab.R2(1) tab.R2(2) ];
    

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end




% check
disp(stat.grip);
groupstat = varfun(@nanmean,stat.grip,'GroupingVariables','group');
disp(groupstat);


%%% 3.2 display
% parameters definition
fig = figure; set(fig,'Name','grip_ftd');
hold on;clear h;
yvar = {'kr','ke','kf','tau','R2_force','R2_velocity'};
nvar = numel(yvar);
xpos = [1:nvar];
xticktext = {'kr','ke','kf','tau','R2_force','R2_velocity'};
ytext = 'parameters (au.)';

% plot
ip=[1:4];
% subplot(1,3,1);hold on;
for ivar=ip
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
ylabel(ytext);

% format
setFigProper('FontSize',20,'LineWidth',2);


%% 9/ learning analysis
%%%% 8.1 model-free behavioral results
%%% .1 statistics
% variable definition
ngroup = 2;ndim=3;
taskname = 'learning';
nqtle = 30;
nbin = [ nqtle,ndim,nsub ];
GROUP = nan(nsub,1);
X = nan(nbin);
Y = nan(nbin);
Yhat = nan(nbin);

% first-level statistics
for isub = 1:nsub
    try
    % select
    tab = data{isub}.battery.table; 
    selection =  (tab.task==taskname) ;
    tab = tab(selection,:);

    % variables
    dim = tab.pairValence;
    nt = tab.trialNumberByValence;
    correct = tab.isOptimalChoice;
    correct2 = tab.predicted_isOptimalChoice;

    % stats
    igroup = (isub>n_control) + 1;
    GROUP(isub)=igroup;
    % 1-select only version with nt=30
    xbin = nt;
    % 2-select all versions 
%     xbin = quantileranks(nt,nqtle);
    X(:,:,isub) = tools.tapply(nt,{xbin,dim},@nanmean);
    Y(:,:,isub)  = tools.tapply(correct,{xbin,dim},@nanmean);
    Yhat(:,:,isub)  = tools.tapply(correct2,{xbin,dim},@nanmean);

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end

% check
% disp(Y);
% disp(Y2);
% disp(Y3);


%%% .2 display
% parameters definition
fig = figure; set(fig,'Name','learning_ftd');
hold on;clear h;


%%%% gain plot
% xticktext = {};
xtext = ' trial number (n)';
ytext = 'choice = correct (%)';
titletext = 'gain';
xlimits = [0 nqtle+1];
ylimits = [0 1];
dsub=3;
subplot(1,2,1);hold on;
% plot
for ivar=1:nvar
    for ig = [1 2]
        % data2plot
        xx = nanmean(X(:,3,GROUP==ig),dsub);
        yy = nanmean(Y(:,3,GROUP==ig),dsub);
        yy2 = nanmean(Yhat(:,3,GROUP==ig),dsub);
        zz = sem(Y(:,3,GROUP==ig),dsub);
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
ax.XTick = [0:3:nqtle] ;
% ax.XTickLabel = {'0.01','0.2','0.50','1','5','20'} ;
ax.XLim = xlimits;       
ax.YLim = ylimits;    
xlabel(xtext); 
ylabel(ytext); 
title(titletext);

%%%% loss plot
% xticktext = {};
xtext = ' trial number (n)';
ytext = 'choice = correct (%)';
titletext = 'loss';
xlimits = [0 nqtle+1];
ylimits = [0 1];
dsub=3;
subplot(1,2,2);hold on;
% plot
for ivar=1:nvar
    for ig = [1 2]
        % data2plot
        xx = nanmean(X(:,1,GROUP==ig),dsub);
        yy = nanmean(Y(:,1,GROUP==ig),dsub);
        yy2 = nanmean(Yhat(:,1,GROUP==ig),dsub);
        zz = sem(Y(:,1,GROUP==ig),dsub);
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
ax.XTick = [0:3:nqtle] ;
% ax.XTickLabel = {'0.01','0.2','0.50','1','5','20'} ;
ax.XLim = xlimits;       
ax.YLim = ylimits;    
xlabel(xtext); 
ylabel(ytext); 
title(titletext);

%% 9.2 model-based analysis
%%% 3.1 statistics
% variable definition
ndim = 2; ngroup = 2;
nbin = [ ndim ,2, ngroup , nsub ];
taskname = 'learning';
Y = nan(nbin);

varnames = {'group','alpha','kr','kp','bm','bp','BCA'};
nvar = numel(varnames);
stat.(taskname) = array2table(nan(nsub,nvar)) ;
stat.(taskname).Properties.VariableNames = varnames;

% first-level statistics
for isub = 1:nsub
    try
    % select
    tab = result{isub}.learning.inferential;

    % stats
    igroup = (isub>n_control) + 1;
    stat.(taskname).group(isub) = igroup;
    
    stat.(taskname){isub,[2:7]} = tab{1,[9:11 13:15]};
    

    % exceptional error
    catch
        fprintf(['statistics not available for sub %d \n'],isub);
    end
end




% check
disp(stat.learning);
groupstat = varfun(@nanmean,stat.learning,'GroupingVariables','group');
disp(groupstat);


%%% 3.2 display
% parameters definition
fig = figure; set(fig,'Name','learning_ftd');
hold on;clear h;
yvar = {'alpha','kr','kp','bm','bp','BCA'};
nvar = numel(yvar);
xpos = [1:nvar];
xticktext = {'alpha','kr','kp','bm','bp','BCA'};
ytext = 'parameters (au.)';

% plot
ip=1;
subplot(1,3,1);hold on;
for ivar=ip
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
% ax.XLim = [0 xpos(end)+1];       
% ax.YLim = [-4 6];       
ylabel(ytext); 
        
% plot
ip=[2 3 4];
subplot(1,3,2);hold on;
for ivar=ip
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

% plot
ip=6;
subplot(1,3,3);hold on;
for ivar=ip
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

% format
setFigProper('FontSize',20,'LineWidth',2);


% format
setFigProper('FontSize',20,'LineWidth',2);

%% choice 1D task
    ndim = 4; ngroup = 2;
    nbin = [ ndim , ngroup , nsub ];

    % variables
        Y = nan(nbin);

        for isub = 1:nsub

            % select
                tab = data{isub}.battery.table; 
                selection =  (tab.task=='choice') ;
                tab = tab(selection,:);

            % variables
                dv = tab.differenceRatingValue/100 ;
                dim = tab.dimension;
                choice = tab.sideChoice;
                correct = double(sign(dv)==sign(choice));
                correct(dv==0) = NaN;

            % stats
                [~,subdim] = ismember(unique(dim),dimensionList);
                igroup = (isub>n_control) + 1;
                ysub = tools.tapply(correct,{dim},@nanmean);
                Y(subdim,igroup,isub) =  ysub;

        end
    
    % display
      fig = figure; set(fig,'Name','choiec_ftd');
            
        dsub=3;
%         y = nanmedian(Y,dsub);
        y = nanmean(Y,dsub);
        z = sem(Y,dsub);
%         z = nanstd(Y,0,dsub);


        hold on
        for ig = [1 2]
            x = [1 4];
            xx = x + (ig-1)*0.6;
            [ h(ig)  ] = barplot( xx ,y(x,ig),z(x,ig), col{ig} );
            h(ig).BarWidth = 0.2;
        end
        legend([h(1) h(2)],cellstr(groupList));

    % legending
        ax = gca; 
        ax.TickLength = [0 0];
        ax.XTick = x +0.5 ;
        ax.XTickLabel = cellstr(dimensionList(x));
%         ax.XLim = [0 x(end)+1];
        ax.YLim = [0.5 1];       
        ylabel('predicted choice (%)'); 
%         ylabel('extreme rating (%)'); 

        
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
%% choice 1D task
    ndim = 3; ngroup = 2;
    nbin = [ ndim , ngroup , nsub ];

    % variables
        Y = nan(nbin);

        for isub = 1:nsub

            % select
                tab = data{isub}.battery.table; 
                selection =  (tab.task=='weight') ;
                tab = tab(selection,:);

            % variables
                dim = tab.dimension;
                accept = tab.isGoChoice;

            % stats
                [~,subdim] = ismember(unique(dim),dimensionList);
                idim = subdim-7;
                igroup = (isub>n_control) + 1;
                ysub = tools.tapply(accept,{dim},@nanmean);
                Y(idim,igroup,isub) =  ysub;

        end
    
    % display
      fig = figure; set(fig,'Name','weight_ftd');
            
        dsub=3;
%         y = nanmedian(Y,dsub);
        y = nanmean(Y,dsub);
        z = sem(Y,dsub);
%         z = nanstd(Y,0,dsub);


        hold on;
        for ig = [1 2]
            x = [1];
            xx = x + (ig-1.5)*0.3;
            [ h(ig)  ] = barplot( xx ,y(x,ig),z(x,ig), col{ig} );
            h(ig).BarWidth = 0.2;
        end
        legend([h(1) h(2)],cellstr(groupList));

    % legending
        ax = gca; 
        ax.TickLength = [0 0];
        ax.XTick = x  ;
        ax.XTickLabel = cellstr(dimensionList(x+7));
        ax.XLim = [x(1) x(end)] + 0.75*[-1 1];
        ax.YLim = [0 1];       
        ylabel('acceptance rate (%)'); 
%         ylabel('extreme rating (%)'); 

        
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);

%% choice IT task
    ndim = 4; ngroup = 2;
    nbin = [ ndim , ngroup , nsub ];

    % variables
        Y = nan(nbin);

        for isub = 1:nsub
            try
            % select
                tab = data{isub}.battery.table; 
                selection =  (tab.task=='discount') ;
                tab = tab(selection,:);

            % variables
                dim = tab.dimension;
                patient = tab.isDelayedItemChoice;

            % stats
                [~,subdim] = ismember(unique(dim),dimensionList);
                idim = subdim;
                igroup = (isub>n_control) + 1;
                ysub = tools.tapply(patient,{dim},@nanmean);
                Y(idim,igroup,isub) =  ysub;
            end
        end
    
    % display
      fig = figure; set(fig,'Name','discount_ftd');
            
        dsub=3;
%         y = nanmedian(Y,dsub);
        y = nanmean(Y,dsub);
        z = sem(Y,dsub);
%         z = nanstd(Y,0,dsub);


        hold on;
        for ig = [1 2]
            x = [1];
            xx = x + (ig-1.5)*0.3;
            [ h(ig)  ] = barplot( xx ,y(x,ig),z(x,ig), col{ig} );
            h(ig).BarWidth = 0.2;
        end
        legend([h(1) h(2)],cellstr(groupList));

    % legending
        ax = gca; 
        ax.TickLength = [0 0];
        ax.XTick = x  ;
        ax.XTickLabel = cellstr(dimensionList(x));
        ax.XLim = [x(1) x(end)] + 0.75*[-1 1];
        ax.YLim = [0 1];       
        ylabel('patient choice (%)'); 
%         ylabel('extreme rating (%)'); 

        
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
    
%% learning task
    ndim = 3; ngroup = 2;
    nbin = [ ndim , ngroup , nsub ];

    % variables
        Y = nan(nbin);

        for isub = 1:nsub
            try
            % select
                tab = data{isub}.battery.table; 
                selection =  (tab.task=='learning') ;
                tab = tab(selection,:);

            % variables
                dim = tab.pairValence;
                correct = tab.isOptimalChoice;

            % stats
                igroup = (isub>n_control) + 1;
                ysub = tools.tapply(correct,{dim},@nanmean);
                Y(:,igroup,isub) =  ysub;
            end
        end
    
    % display
      fig = figure; set(fig,'Name','learning_ftd');
            
        dsub=3;
%         y = nanmedian(Y,dsub);
        y = nanmean(Y,dsub);
        z = sem(Y,dsub);
%         z = nanstd(Y,0,dsub);


        hold on;
        for ig = [1 2]
            x = [1 3];
            xx = x + (ig-1.5)*0.3;
            [ h(ig)  ] = barplot( xx ,y(x,ig),z(x,ig), col{ig} );
            h(ig).BarWidth = 0.2;
        end
        chance = [0.5 0.5];
        plot([0 xx(end)+2],chance,'k--');
        legend([h(1) h(2)],cellstr(groupList));

    % legending
        xNames = {'Loss','Gain'};
    
        ax = gca; 
        ax.TickLength = [0 0];
        ax.XTick = x  ;
        ax.XTickLabel = xNames;
        ax.XLim = [x(1) x(end)] + 0.75*[-1 1];
        ax.YLim = [0 1];       
        ylabel('correct choice (%)'); 

        
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);

%% learning task parameters
%     nparam = 5;
    nparam = 3; 
    ngroup = 2;
    nbin = [ nparam , ngroup , nsub ];

    % variables
        Y = nan(nbin);

        for isub = 1:nsub
            try
            % select
                tab = result{isub}.learning.inferential  ; 

            % variables
                paramNames = {'alpha','kR','kP','bm','bp'};
                paramNames = {'repetition_previousLoss','repetition_Neutral','repetition_previousGain'};

                param = [];
                for ivar = 1:numel(paramNames)
                   param = [param tab.(paramNames{ivar}) ] ;
                end

            % stats
                igroup = (isub>n_control) + 1;
                Y(:,igroup,isub) =  param;
            end
        end
    
    % display
      fig = figure; set(fig,'Name','learning_param_ftd');
            
        dsub=3;
%         y = nanmedian(Y,dsub);
        y = nanmean(Y,dsub);
        z = sem(Y,dsub);
%         z = nanstd(Y,0,dsub);


        hold on;
        for ig = [1 2]
            x = [1:nbin(1)];
%             x = [1 2 3];
            xx = x + (ig-1.5)*0.3;
            [ h(ig)  ] = barplot( xx ,y(x,ig),z(x,ig), col{ig} );
            h(ig).BarWidth = 0.2;
%             [ h(ig) ] = myboxplot( xx , reshape(Y(x,ig,:),numel(x),nsub)' , col{ig} );
            
        end
        legend([h(1) h(2)],cellstr(groupList));

    % legending
        paramNames = {'repetition|loss_{t-1}','repetition|neutral_{t-1}','repetition|gain_{t-1}'};
        xNames = paramNames;
    
        ax = gca; 
        ax.TickLength = [0 0];
        ax.XTick = x  ;
        ax.XTickLabel = xNames;
        ax.XLim = [x(1) x(end)] + 0.75*[-1 1];
%         ax.YLim = [0 1];       

        
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
%% effort task parameters
    nparam = 4; ngroup = 2;
    nbin = [ nparam , ngroup , nsub ];

    % variables
        Y = nan(nbin);

        for isub = 1:nsub
            try
            % select
                tab = result{isub}.grip.inferential  ; 

            % variables
                paramNames = {'kR','kE','kF','fmax'};
                param = [];
                for ivar = 1:numel(paramNames)
                   param = [param tab.(paramNames{ivar}) ] ;
                end

            % stats
                igroup = (isub>n_control) + 1;
                Y(:,igroup,isub) =  param;
            end
        end
    
    % display
      fig = figure; set(fig,'Name','effort_param_ftd');
            
        dsub=3;
%         y = nanmedian(Y,dsub);
        y = nanmean(Y,dsub);
        z = sem(Y,dsub);
%         z = nanstd(Y,0,dsub);


        hold on;
        for ig = [1 2]
            x = [1:nbin(1)];
            x = [1 2 3];
            xx = x + (ig-1.5)*0.3;
%             [ h(ig)  ] = barplot( xx ,y(x,ig),z(x,ig), col{ig} );
%             h(ig).BarWidth = 0.2;
            [ h(ig) ] = myboxplot( xx , reshape(Y(x,ig,:),numel(x),nsub)' , col{ig} );
            
        end
        legend([h(1) h(2)],cellstr(groupList));

    % legending
        xNames = paramNames;
    
        ax = gca; 
        ax.TickLength = [0 0];
        ax.XTick = x  ;
        ax.XTickLabel = xNames;
        ax.XLim = [x(1) x(end)] + 0.75*[-1 1];
%         ax.YLim = [0 1];       

        
    % reformat
    setFigProper('FontSize',20,'LineWidth',2);
    
    
%% heart rate  parameters
    
