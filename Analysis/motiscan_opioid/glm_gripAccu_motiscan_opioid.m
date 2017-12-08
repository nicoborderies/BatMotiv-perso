function [p,coefNames,stat] = glm_gripAccu_motiscan_opioid(data,yname,varnames,formula,varargin)

% check options
optionList = {'fname','ftransform','sessionSplit','order','statisticalTest','plotType'};
defaultList = { @nanmean, @identity,0,[],{'drug_ttest'},{'bar'}};
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
for i = 1:numel(statisticalTest)
    switch statisticalTest{i}
        case 'drug_ttest'
            col = { [0 0 1]*0.75 , [1 1 1]*0.5 , [1 0 0]*0.75 };
        case 'opioid_ttest'
            col = { [1 1 1]*0.5 , [1 0 1]*0.75 };
    end
end
%%% transient 
taskname = 'gripAccu';
nbin = [ ntrt , nsub ];


% prepare variables
Y = []; % conditional means by treatment
S = nan(nbin);
R2 = [];

% group averaging
for isub = 1:nsub % subject loop
    % select
%         tab = data{isub}.battery.table; 
%         selection =  (tab.task==taskname) ;
%         tab = tab(selection,:);
        tab = data{isub}.(taskname).table; 

    % variables
        y = tab.(yname);
        y = ftransform(y);
        trt = tab.treatment;
        trt = removecats(trt,'0');
        trt = reordercats(trt,treatmentList);
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
        % predictive variables
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
%         disp(coef);
    
        [~,subtrt] = ismember(unique(trt),treatmentList);
        sess = splitapply(fname,session,double(trt));
    % store
        Y =  [Y,coef.Estimate] ;
        R2 = [R2,stat.Rsquared.Ordinary];
        S(subtrt,isub) =  sess;
        if isub==1; coefNames = coef.Properties.RowNames; end
end
%%% average
dsub=2;
ncoef = numel(coefNames);
y = nanmean(Y,dsub);
r2 = nanmean(R2);
%%% concatenate
Y2 = reshape(Y,nsub*ncoef,1);
COEF = repmat(coefNames,nsub,1);
COEF = nominal(COEF);
COEF = reordercats(COEF,coefNames);

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
for i = 1:numel(statisticalTest)
    switch statisticalTest{i}
        case 'drug_ttest'
            drugLabel = reordercats(drugLabel,treatmentList);
        case 'opioid_ttest'
            treatmentList = {'placebo','opioid'};
            drugLabel = reordercats(drugLabel,treatmentList);
    end
end
COEFLABEL = repmat(coefLabel,nsub,1);
DRUGLABEL = repmat(drugLabel,nsub,1);

% statistics 2nd level
%%% ttest
[h,p,~,stat] = ttest(Y');


% display
clear g;
alpha=0.8;
g = gramm('x',COEFLABEL,'y',Y2,'color',DRUGLABEL);
g.set_color_options('map',vertcat(col{:}),'lightness',100);
g.set_order_options('color',treatmentList);
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
g.set_names('x',xname,'y','coefficients');
g.set_title(['GLM : ' yname ]);
if sessionSplit || ~isempty(order); g.axe_property('XLim',[0 4]); end
ylimits = [ min(mean(Y,2)-2*sem(Y,2)) max(mean(Y,2)+2*sem(Y,2))];
g.axe_property('YLim',ylimits);
g.axe_property('XTickLabelRotation',45);
g.draw; 
axes(g.facet_axes_handles);
hold on; plot(gca,[0:ncoef+1],zeros(1,ncoef+2),'k');
xpos = double(nominal(coefLabel)) + (double(nominal(drugLabel))-2)*0.2;
sigstar( num2cell([xpos]),p); % display significance
set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                    'LineWidth',1.5,'FaceAlpha',alpha);

end