function [p,score,stat] = barcomp_stat_motiscan_opioid(result,varnames,varlegend,taskname,varargin)

% check options
optionList = {'fname','ftransform','sessionSplit','order','statisticalTest','plotType'};
defaultList = { @nanmean, @identity,0,[],...
                {'naloxone_contrast','morphine_contrast'},...
                {'bar'}};
[fname,ftransform,sessionSplit,order,statisticalTest,plotType] = check_options(varargin,optionList,defaultList);


% parameters
%%% lists
dimensionList = nominal({'writtenReward','visualReward','writtenPunishment','writtenEffort',...
                            'monetaryReward','monetaryPunishment','gripEffort',...
                            'RewardEffort','PunishmentEffort','RewardPunishment',...
                            'Gain','Loss','GainLoss','GainEmotion'});
taskList = nominal( {'rating','choice','weight','discount',...
                    'grip','gripIAPS','gripAccu','mental','learning'});
treatmentList = {'naloxone','placebo','morphine'};
ntrt = numel(treatmentList);
nsub = numel(result);
%%% display
col = { [0 0 1]*0.75 , [1 1 1]*0.5 , [1 0 0]*0.75 };
%%% transient 
nbin = [ ntrt , nsub ];


% prepare variables
Y = []; % conditional means by treatment
Y2 = []; % mean difference to placebo
S = [];
T = [];
V = [];

% group averaging
for isub = 1:nsub % subject loop
    % select
        tab = result{isub}.(taskname).inferential;
    % variables
        y = [];
        for j = 1:numel(varnames)
            y = [y, tab.(varnames{j})];
        end
        y = ftransform(y);
        Y =  [Y ; y];
        T = [T;tab.treatment];
        S = [S;[1 2 3]'];
        V = [V; repmat(nominal(varnames),3,1)];

    % normalization
        y2 = y - nanmean(y);
%         y2 = y - y(tab.treatment=='placebo',:);
        Y2 =  [Y2 ; y2];

end

%%% average
Y2 = Y2 + nanmean(Y);
% Y2 = Y2 + nanmean(Y(T=='placebo',:));

%%% concatenate
if ~isempty(order)
    O = repmat(order',3,1);  
    O = reshape(O,nsub*ntrt,1);
else
    O = repmat(1,nsub*ntrt,1);
end


% statistics 2nd level
for i = 1:numel(varnames)
    %%% contrast
    predictor = [ (T=='naloxone') , (T=='placebo'), (T=='morphine') , S, O ];
    y = Y2(:,i); 
    if sessionSplit
        formula = [varnames{i} ' ~ -1 + placebo + naloxone + morphine + session + session:morphine + session:naloxone'];
        stat{i} = fitglm(predictor,y,...
                       formula,'VarNames',{'naloxone','placebo','morphine','session','order',varnames{i}});
        coef = stat{i}.Coefficients;
        disp(coef);
    elseif ~isempty(order)
        formula = [varnames{i} ' ~ -1 + placebo + naloxone + morphine + session + order + order:morphine + order:naloxone'];
        stat{i} = fitglm(predictor,y,...
                       formula,'VarNames',{'naloxone','placebo','morphine','session','order',varnames{i}});
        coef = stat{i}.Coefficients;
        disp(coef);
    else
        switch statisticalTest{1}
            case 'opioid_ttest'
                predictor = [ (T=='naloxone')*(-1) + (T=='placebo')*0 + (T=='morphine')*1 , S, O ];
                formula = [varnames{i} ' ~ 1 + opioid + session '];
                stat{i} = fitglm(predictor,y,...
                            formula,'VarNames',{'opioid','session','order',varnames{i} });
                coef = stat{i}.Coefficients;
                disp(coef);
            otherwise
                formula = [varnames{i} ' ~ -1 + placebo + naloxone + morphine + session '];
                stat{i} = fitglm(predictor,y,...
                               formula,'VarNames',{'naloxone','placebo','morphine','session','order',varnames{i}});
                coef = stat{i}.Coefficients;
                disp(coef);
        end
    end
end

%%% statistical test
p = nan(numel(statisticalTest),numel(varnames));
score = nan(numel(statisticalTest),numel(varnames));
for j = 1:numel(varnames)
    for i = 1:numel(statisticalTest)
        switch statisticalTest{i}
            case 'morphine_contrast'
                contrast = [0 -1 1];
            case 'naloxone_contrast'
                contrast = [1 -1 0];
            case 'opioid_contrast'
                contrast = [0 -1 1; -1 1 0; -1 0 1];
            case 'opioid_ttest'
                contrast = [0 1];
            case 'naloxone_ttest'
                dy = Y2(T=='naloxone')-Y2(T=='placebo');
                if mean(dy)>0; side = 'right';
                else side = 'left';
                end
                [~,p(i,j),~,stat] = ttest(dy,0,'Tail',side);
                score(i,j) = stat.tstat;
            case 'morphine_ttest'
                dy = Y2(T=='morphine')-Y2(T=='placebo');
                if  mean(dy)>0; side = 'right';
                else side = 'left';
                end
                [~,p(i,j),~,stat] = ttest(dy,0,'Tail',side);
                score(i,j) = stat.tstat;
        end

        if sessionSplit
            ncofactor = 3;
        elseif ~isempty(order)
            ncofactor = 4;
        elseif strcmp(statisticalTest{i},'opioid_ttest')
            ncofactor = 1;
        else
            ncofactor = 1;
        end
%         try
            contrast = [contrast , zeros(size(contrast,1),ncofactor) ];
            [p(i,j),score(i,j),d] = coefTest(stat{j},contrast)
%         end
    end
end

% display
%%% concatenate
Y2 = reshape(Y2,nsub*ntrt*numel(varnames),1);
V = reshape(V,nsub*ntrt*numel(varnames),1);
T = repmat(T,1,numel(varnames));
T = reshape(T,nsub*ntrt*numel(varnames),1);

clear g;
alpha=0.8;
if sessionSplit
    g = gramm('x',S,'y',Y2,'color',T);
elseif ~isempty(order)
    g = gramm('x',O,'y',Y2,'color',T);
else
    g = gramm('x',V,'y',Y2,'color',T);
end
g.set_color_options('map',vertcat(col{:}),'lightness',100);
if numel(unique(V))>1
    g.set_order_options('x',(varnames),'color',treatmentList);
else
    g.set_order_options('color',treatmentList);
end
for i = 1:numel(plotType)
    switch plotType{i}
        case 'dot'
            g.geom_point();
        case 'jitter'
            g.geom_jitter('height',0.01);
        case 'bar'
            g.stat_summary('type','sem','geom',{'bar','black_errorbar'},'width',0.9);
        case 'boxplot'
            g.stat_boxplot('width',0.9);
        case 'violin'
            g.stat_violin();
    end
end
if sessionSplit==1
    xname = 'session number';
elseif ~isempty(order)
    xname = 'task order';
else
    xname = '';
end
if numel(unique(V))>1
    g.set_names('x',xname,'y',varlegend,'color','treatment');
else
    g.set_names('x',xname,'y',varnames{1},'color','treatment');
end
if sessionSplit || ~isempty(order); g.axe_property('XLim',[0 4]); end
g.axe_property('YLim',[mean(min(Y2)) mean(max(Y2))] + [-0.2 0.2]*mean(mean(Y2)));
g.axe_property('XTickLabelRotation',45);
g.draw;
axes(g.facet_axes_handles);
% display significance
for j = 1:numel(varnames)
    for i = 1:numel(statisticalTest)
        switch statisticalTest{i}
            case 'morphine_contrast'
                xloc = [0 1];
            case 'naloxone_contrast'
                xloc = [-1 0];
            case 'opioid_contrast'
                xloc = [0 0];
            case 'opioid_ttest'
                xloc = [0 0];
            case 'naloxone_ttest'
                xloc = [0 1];
            case 'morphine_ttest'
                xloc = [-1 0];
        end
        xloc = j + 0.2*xloc;
        pstar = p(i,j); if pstar>0.10; pstar=NaN; end
        sigstar({xloc},pstar); 
    end
end
set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                    'LineWidth',1.5,'FaceAlpha',alpha);

end