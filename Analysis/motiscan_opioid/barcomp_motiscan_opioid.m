function [p,score,g] = barcomp_motiscan_opioid(data,varname,varlegend,taskname,varargin)

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
nsub = numel(data);
%%% display
col = { [0 0 1]*0.75 , [1 1 1]*0.5 , [1 0 0]*0.75 };
%%% transient 
nbin = [ ntrt , nsub ];


% prepare variables
Y = nan(nbin); % conditional means by treatment
Y2 = nan(nbin); % mean difference to placebo
S = nan(nbin);

% group averaging
for isub = 1:nsub % subject loop
    % select
        tab = data{isub}.(taskname).table; 
        selection =  (tab.task==taskname) ;
        tab = tab(selection,:);
    % variables
        y = tab.(varname);
        y = ftransform(y);
        trt = tab.treatment;
        trt = removecats(trt,'0');
        trt = reordercats(trt,treatmentList);
        session = tab.sessionNumber;
        
    % stats
        [~,subtrt] = ismember(unique(trt),treatmentList);
        ysub = splitapply(fname,y,double(trt));
        sess = splitapply(fname,session,double(trt));
    % store
        Y(subtrt,isub) =  ysub;
        S(subtrt,isub) =  sess;
end
%%% average
dsub=2;
y = nanmean(Y,dsub);
z = sem(Y2,dsub);
%%% concatenate
% Y2 = Y-Y(2,:)+nanmean(Y(2,:)); % normalize by the subject mean under placebo
Y2 = Y - nanmean(Y,1) + nanmean(nanmean(Y)); % normalize by the subject mean
Y = reshape(Y,nsub*ntrt,1);
Y2 = reshape(Y2,nsub*ntrt,1);
T = repmat(nominal(treatmentList'),nsub,1);
T = reordercats(T,treatmentList);

S = reshape(S,nsub*ntrt,1);
SUB = repmat([1:nsub],ntrt,1);
SUB = reshape(SUB,nsub*ntrt,1);

if ~isempty(order)
    O = repmat(order',3,1);  
    O = reshape(O,nsub*ntrt,1);
else
    O = ones(nsub*ntrt,1);
end


% statistics 2nd level
%%% GLM
predictor = [ (T=='naloxone') , (T=='placebo'), (T=='morphine') , S, O ];
y = Y2; 
if sessionSplit
    formula = [varname ' ~ -1 + placebo + naloxone + morphine + session + session:morphine + session:naloxone'];
    stat = fitglm(predictor,y,...
                   formula,'VarNames',{'naloxone','placebo','morphine','session','order',varname});
    coef = stat.Coefficients;
    disp(coef);
elseif ~isempty(order)
    formula = [varname ' ~ -1 + placebo + naloxone + morphine + session + order + order:morphine + order:naloxone'];
    stat = fitglm(predictor,y,...
                   formula,'VarNames',{'naloxone','placebo','morphine','session','order',varname});
    coef = stat.Coefficients;
    disp(coef);
else
    formula = [varname ' ~ -1 + placebo + naloxone + morphine + session '];
    stat = fitglm(predictor,y,...
                   formula,'VarNames',{'naloxone','placebo','morphine','session','order',varname});
    coef = stat.Coefficients;
    disp(coef);
end

%%% statistical test
p = nan(1,numel(statisticalTest));
score = nan(1,numel(statisticalTest));

for i = 1:numel(statisticalTest)
    switch statisticalTest{i}
        case 'morphine_contrast'
            contrast = [0 -1 1];
        case 'naloxone_contrast'
            contrast = [1 -1 0];
        case 'opioid_contrast'
%             contrast = [0 -1 1; -1 1 0; -1 0 1];
            contrast = [-1 0 1];
        case 'opioid_ttest'
            predictor = [ (T=='naloxone')*(-1) + (T=='placebo')*0 + (T=='morphine')*1 , S, O ];
            formula = [varname ' ~ 1 + opioid + session '];
            stat = fitglm(predictor,y,...
                        formula,'VarNames',{'opioid','session','order',varname});
            coef = stat.Coefficients;
            disp(coef);
            contrast = [0 1];
        case 'naloxone_ttest'
            dy = Y2(T=='naloxone')-Y2(T=='placebo');
            if mean(dy)>0; side = 'right';
            else side = 'left';
            end
            [~,p(i),~,stat] = ttest(dy,0,'Tail',side);
            score(i) = stat.tstat;
        case 'morphine_ttest'
            dy = Y2(T=='morphine')-Y2(T=='placebo');
            if  mean(dy)>0; side = 'right';
            else side = 'left';
            end
            [~,p(i),~,stat] = ttest(dy,0,'Tail',side);
            score(i) = stat.tstat;
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
    try
        contrast = [contrast , zeros(size(contrast,1),ncofactor) ];
        [p(i),score(i),d] = coefTest(stat,contrast)
    end
end


% display
clear g;
alpha=0.7;   
if sessionSplit
    g = gramm('x',S,'y',Y2,'color',T);
elseif ~isempty(order)
    g = gramm('x',O,'y',Y2,'color',T);
else
    g = gramm('x',T,'y',Y2,'color',T);
end
g.set_color_options('map',vertcat(col{:}),'lightness',100);
g.set_order_options('x',treatmentList,'color',treatmentList);
for i = 1:numel(plotType)
    switch plotType{i}
        case 'dot'
            g.geom_point();
        case 'jitter'
            g.geom_jitter('height',0.01);
        case 'line'
            for isub=1:nsub
%                 indsub = [1:ntrt]' + (isub-1)*ntrt ; 
                g.update('subset',(SUB==isub),'color',[]);
                g.set_color_options('map',[0 0 0]);
                g.geom_line();
%                 g.geom_abline('style','k-');
                g.draw();
                ysub = splitapply(@nanmean,Y(SUB==isub),double(T(SUB==isub)));
                if ysub(3)-ysub(1)>=0; lstyle = '-';else;lstyle = '--';end
                set(g.results.geom_line_handle,'Color',0*[1 1 1],'LineStyle',lstyle,'LineWidth',1)
            end
            g.update('x',T,'y',Y2,'color',T,'subset',[]);
            g.set_color_options('map',vertcat(col{:}),'lightness',100);
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
g.set_names('x',xname,'y',varlegend,'color','treatment');
if sessionSplit || ~isempty(order); g.axe_property('XLim',[0 4]); end
g.axe_property('YLim',[min(Y2) max(Y2)] + [-0.2 0.2]*mean(Y2));
g.draw;
axes(g.facet_axes_handles);
% display significance
for i = 1:numel(statisticalTest)
    switch statisticalTest{i}
        case 'morphine_contrast'
            xloc = [2 3];
        case 'naloxone_contrast'
            xloc = [1 2];
        case 'opioid_contrast'
            xloc = [2 2];
        case 'opioid_ttest'
            xloc = [2 2];
        case 'naloxone_ttest'
            xloc = [1 2];
        case 'morphine_ttest'
            xloc = [2 3];
    end
    pstar = p(i); if pstar>0.10; pstar=NaN; end
    sigstar({xloc},pstar); 
end
set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                    'LineWidth',1.5,'FaceAlpha',alpha);

end