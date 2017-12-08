%% motiscan_opioid_2level_gripAccu

%% 1. Effect of drugs on major measurements

% define test from dataset
taskname = 'gripAccu';

%% parameters
% define metrics
varname = 'forceDuration';
ftransform = @ (x) normalize(x,'zscore');
% ftransform = @identity;
fname = @nanmean;
varlegend = 'force duration, drug-placebo (sec)';

% define statistical procedure
statisticalTest = {'opioid_contrast'};

% display options
% - plot type
plotType = {'jitter','boxplot'};

%% aggregate data
% - prepare data
nbin = [ ntrt , nsub ];
Y = nan(nbin); % conditional means by treatment
Y2 = nan(nbin); % mean difference to placebo
S = nan(nbin);
% - condition data
for isub = 1:nsub % subject loop
    % select
        tab = data{isub}.(taskname).table; 
        selection =  (tab.task==taskname) ;
        tab = tab(selection,:);
    % variables
        y = tab.(varname);
        y2 = ftransform(y);
        trt = tab.treatment;
        trt = removecats(trt,'0');
        trt = reordercats(trt,treatmentList);
        session = tab.sessionNumber;
        
    % stats
        [~,subtrt] = ismember(unique(trt),treatmentList);
        ysub = splitapply(fname,y,double(trt));
        ysub2 = splitapply(fname,y2,double(trt));
        sess = splitapply(fname,session,double(trt));
    % store
        Y(subtrt,isub) =  ysub;
        Y2(subtrt,isub) =  ysub2;
        S(subtrt,isub) =  sess;
end
% - concatenate data
% Y2 = Y2 - Y2(2,:); % normalize by the subject mean under placebo
Y2 = Y2 - nanmean(Y2,1); % normalize by the subject mean
Y = reshape(Y,nsub*ntrt,1);
Y2 = reshape(Y2,nsub*ntrt,1);
T = repmat(nominal(treatmentList'),nsub,1);
T = reordercats(T,treatmentList);
S = reshape(S,nsub*ntrt,1);
SUB = repmat([1:nsub],ntrt,1);
SUB = reshape(SUB,nsub*ntrt,1);
O = ones(nsub*ntrt,1);

%% 2level statistics
% construct linear model
predictor = [ (T=='naloxone') , (T=='placebo'), (T=='morphine') , S, O ];
% predictor = predictor(T~='placebo',:);
y = Y2; 
% y = y(T~='placebo'); 
formula = [varname ' ~ -1 + placebo + naloxone + morphine + session '];
stat = fitglm(predictor,y,...
               formula,'VarNames',{'naloxone','placebo','morphine','session','order',varname});
coef = stat.Coefficients;
disp(coef);
% statistical inference
p = nan(1,numel(statisticalTest));
score = nan(1,numel(statisticalTest));
contrast = [-1 0 1];
for i = 1:numel(statisticalTest)
%     try
        ncofactor = 1;
        contrast = [contrast , zeros(size(contrast,1),ncofactor) ];
        [p(i),score(i),d] = coefTest(stat,contrast)
        
%         dy = Y2(T=='morphine')-Y2(T=='naloxone');
%         [~,p(i),~,stat] = ttest(dy,0)
%         score(i) = stat.tstat
        
%     end

end

%% display
fig = figure;
clear g;
alpha=0.7;  
g(1,1) = gramm('x',T,'y',Y,'color',T,'subset',(T=='placebo'));
g(1,1).set_color_options('map',vertcat(col{[2]}),'lightness',100);
g(1,1).set_order_options('x',treatmentList,'color',treatmentList);
g(1,1).geom_jitter('height',0.01);
g(1,1).stat_boxplot('width',0.3);
g(1,1).set_names('x','','y','force duration (sec)','color','treatment');
g(1,1).axe_property('YLim',[min(Y) max(Y)] + [-0.2 0.2]*mean(Y));
g(1,1).axe_property('XLim',[0.5 1.5]);

g(1,2) = gramm('x',T,'y',Y2,'color',T,'subset',(T~='placebo'));
g(1,2).set_color_options('map',vertcat(col{[1 3]}),'lightness',100);
g(1,2).set_order_options('x',treatmentList([1 3]),'color',treatmentList([1 3]));
g(1,2).geom_jitter('height',0.01);
g(1,2).stat_boxplot('width',0.9);
g(1,2).set_names('x','','y','force duration, drug-average (z-score)','color','treatment');
g(1,2).axe_property('YLim',[min(Y2) max(Y2)] + [-0.2 0.2]*mean(Y2));

g.draw;
axes(g(1,2).facet_axes_handles);
for i = 1:numel(statisticalTest)
    xloc = [1 2];
    pstar = p(i); if pstar>0.10; pstar=NaN; end
    sigstar({xloc},pstar); 
end
set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                    'LineWidth',1.5,'FaceAlpha',alpha,'Interpreter','tex');
l = plot(g(1,2).facet_axes_handles.XLim,[0 0],'--k');
fig.Position = [ 100 100 1000 400];

%% 2. Effect of drugs interaction with experimental factors

%% parameters
% define metrics
varname = 'forceDuration';
ftransform = @ (x) normalize(x,'zscore');
% ftransform = @identity;
fname = @nanmean;
varlegend = 'force duration, drug-placebo (sec)';

% define statistical procedure
statisticalTest = {'drug_ttest'};
%   - model equation
formula = [ varname ' ~ 1 + incentive + difficulty + instruction:difficulty + nblock + ntrial_block + hand',...
           '+ naloxone:(1 + incentive + difficulty + instruction:difficulty + nblock + ntrial_block + hand)',...
           '+ morphine:(1 + incentive + difficulty + instruction:difficulty + nblock + ntrial_block + hand)',...
           '+ session'];
varnames = {'incentive','difficulty','instruction','effort_t',...
            'nblock','ntrial_block','hand',...
            'naloxone','placebo','morphine','opioid','session',varname};
% display options
% - plot type
plotType = {'jitter','boxplot'};

%% aggregate data
% - prepare data
nbin = [ ntrt , nsub ];
Y = []; % conditional means by treatment
S = nan(nbin);
R2 = [];
% - condition data
for isub = 1:nsub % subject loop
    % select
        tab = data{isub}.(taskname).table; 
        selection =  (tab.task==taskname) ;
        tab = tab(selection,:);
    % variables
        y = tab.(varname);
        y2 = ftransform(y);
        trt = tab.treatment;
        trt = removecats(trt,'0');
        trt = reordercats(trt,treatmentList);
        session = tab.sessionNumber;
        naloxone = (trt=='naloxone');
        placebo = (trt=='placebo');
        morphine = (trt=='morphine');
        opioid = (-1)*naloxone + 0*placebo + 1*morphine;
        session = tab.sessionNumber;
        effort = tab.effortDuration;
        rest = tab.restDuration;
        difficulty = tab.costValue;
        incentive = tab.incentiveValue;
        instruction = tab.explicitCost;
        nblock = ceil(tab.trialNumber/8);
        ntrial_block = mod(tab.trialNumber-1,8)+1;
        hand = tab.handSide;
        effort_t = nan(size(y));
        effort_t(2:end) = effort(1:end-1);
        effort_t(ntrial_block<=1) = NaN;
        
    % stats
        predictor = [];
        for i = 1:(numel(varnames)-1)
            eval([ 'x = ' varnames{i} ';' ]);
            if ~islogical(x)
               x = nanzscore(x) ;
            end
            predictor = [predictor,x];
        end
        % fit
        stat = fitglm(predictor,y,...
                       formula,'VarNames',varnames);
        coef = stat.Coefficients;
        [~,subtrt] = ismember(unique(trt),treatmentList);
        sess = splitapply(fname,session,double(trt));
        
    % store
        Y =  [Y,coef.Estimate] ;
        R2 = [R2,stat.Rsquared.Ordinary];
        S(subtrt,isub) =  sess;
        if isub==1; coefNames = coef.Properties.RowNames; end
end
% - concatenate data
dsub=2;
ncoef = numel(coefNames);
y = nanmean(Y,dsub);
r2 = nanmean(R2);
%%% concatenate
Y2 = reshape(Y,nsub*ncoef,1);
COEF = repmat(coefNames,nsub,1);
COEF = nominal(COEF);
COEF = reordercats(COEF,coefNames);
% -- name formating
drugLabel = repmat(nominal('placebo'),numel(coefNames),1);
coefLabel = coefNames;
for i = 1:numel(statisticalTest)
    switch statisticalTest{i}
        case 'drug_ttest'
            drugNames = {'morphine','naloxone'};
        case 'opioid_ttest'
            drugNames = {'opioid'};
    end
end
for i=1:numel(drugNames)
    drugname = drugNames{i};
    finddrug = @(field) contains(field,drugname);
    drugflag = cellfun(finddrug,coefNames);
    drugLabel(drugflag) = nominal(drugname);
    replabel = @(field) strrep(field,[':' drugname ],'');
    coefLabel = cellfun(replabel,coefLabel,'UniformOutput',0);
    replabel = @(field) strrep(field,[drugname],'(Intercept)');
    coefLabel = cellfun(replabel,coefLabel,'UniformOutput',0);
end
COEFLABEL = nominal(repmat(coefLabel,nsub,1));
coefLabel = reordercats(nominal(coefLabel),...
{'(Intercept)','incentive','difficulty','difficulty:instruction','ntrial_block','nblock','hand','session'});
DRUGLABEL = nominal(repmat(drugLabel,nsub,1));
drugLabel = reordercats(nominal(drugLabel),treatmentList);


%% 2level statistics
[h,p,~,stat] = ttest(Y');

%% display
fig = figure;
clear g;
alpha=0.8;
g(1,1) = gramm('x',COEFLABEL,'y',Y2,'color',DRUGLABEL,'subset',DRUGLABEL=='placebo');
g(1,1).set_color_options('map',vertcat(col{2}),'lightness',100);
g(1,1).set_order_options('color',treatmentList);
% g(1,1).geom_jitter('height',0.01);
g(1,1).stat_boxplot('width',0.9);
g(1,1).set_names('x','','y','coefficients');
g(1,1).set_title(['GLM : ' yname ]);
ylimits = [ min(mean(Y,2)-2*sem(Y,2)) max(mean(Y,2)+2*sem(Y,2))];
g(1,1).axe_property('YLim',ylimits);
g(1,1).axe_property('XTickLabelRotation',45);


g.draw; 
axes(g(1,1).facet_axes_handles);
hold on; plot(gca,[0:ncoef+1],zeros(1,ncoef+2),'k--');
xpos = double(nominal(coefLabel)) + (double(nominal(drugLabel))-2)*0.2;
sigstar( num2cell([xpos]),p); % display significance
set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                    'LineWidth',1.5,'FaceAlpha',alpha);
fig.Position = [ 100 100 1200 500];


