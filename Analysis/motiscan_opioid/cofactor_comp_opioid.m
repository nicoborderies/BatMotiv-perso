function [] = cofactor_comp_opioid(data,factor,yname,legend,taskname,varargin)

% check options
optionList = {'fname','ftransform','plotType'};
defaultList = { @nanmean, @identity,...
                {'jitter'}};
[fname,ftransform,plotType] = check_options(varargin,optionList,defaultList);


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


% group averaging
nbin = [ ntrt , nsub ];
Y = nan(nbin); % conditional means by treatment
Y2 = nan(nbin); % mean difference to placebo
S = nan(nbin);
for isub = 1:nsub % subject loop
  % select
        tab = data{isub}.(taskname).table; 
        selection =  (tab.task==taskname) ;
        tab = tab(selection,:);
    % variables
        y = tab.(yname);
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

% normalize
Y2 = Y - nanmean(Y,1) + nanmean(nanmean(Y)); % normalize by the subject mean
% concatenate
Y = reshape(Y,nsub*ntrt,1);
Y2 = reshape(Y2,nsub*ntrt,1);
T = repmat(nominal(treatmentList'),nsub,1);
S = reshape(S,nsub*ntrt,1);
SUB = repmat([1:nsub],ntrt,1);
SUB = reshape(SUB,nsub*ntrt,1);
FACTOR = reshape(repmat(factor,1,ntrt)',nsub*ntrt,1);


% display
clear g;
alpha=0.3;
g = gramm('x',FACTOR,'y',Y2,'color',T);
g.set_color_options('map',vertcat(col{:}),'lightness',100);
g.set_order_options('color',treatmentList);
for i = 1:numel(plotType)
    switch plotType{i}
        case 'dot'
            g.geom_point();
        case 'jitter'
            g.geom_jitter('height',0.01);
        case 'line'
            g.geom_line();
        case 'fit'
            g.stat_glm();
    end
end
g.set_names('x',legend{1},'y',legend{2},'color','treatment');
g.axe_property('YLim',[nanmin(Y2) nanmax(Y2)] + [-0.2 0.2]*nanmean(Y2));
g.draw;
axes(g.facet_axes_handles);
set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                    'LineWidth',1.5,'FaceAlpha',alpha);


end
