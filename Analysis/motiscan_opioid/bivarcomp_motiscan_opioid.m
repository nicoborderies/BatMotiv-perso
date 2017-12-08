function [] = bivarcomp_motiscan_opioid(data,xfactor,yfactor,legend,taskname,varargin)

% check options
optionList = {'fname','ftransform','sessionSplit','order','plotType','interactionfactor'};
defaultList = { @nanmean, @identity,0,[],...
                {'errorbar'},[]};
[fname,ftransform,sessionSplit,order,plotType,interactionfactor] = check_options(varargin,optionList,defaultList);


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
X = [];
Y = []; % conditional means by treatment
Z = [];
Y2 = []; % mean difference to placebo
S = [];
T = [];

% group averaging
for isub = 1:nsub % subject loop
    % select
        tab = data{isub}.(taskname).table; 
%         selection =  (tab.task==taskname) ;
%         tab = tab(selection,:);
    % variables
        y = tab.(yfactor);
        y = ftransform(y);
        x = tab.(xfactor);
        if ~isempty(interactionfactor)
           z =  tab.(interactionfactor);
        end
        trt = tab.treatment;
        trt = removecats(trt,'0');
        trt = reordercats(trt,treatmentList);
        session = tab.sessionNumber;
        
    % stats
        [~,subtrt] = ismember(unique(trt),treatmentList);
        if isempty(interactionfactor)
            [g,gx,gt] = findgroups(x,(trt));
            ysub = splitapply(fname,y,g);
            y2 = ysub;
            for i=unique(gx)'
                y2(gx==i) = y2(gx==i) - mean(y2(gx==i));
            end
        else
            [g,gx,gz,gt] = findgroups(x,z,(trt));
            ysub = splitapply(fname,y,g);
            y2 = ysub;
            [g] = findgroups(gx,gz);
            for i=unique(g)'
                y2(g==i) = y2(g==i) - mean(y2(g==i));
            end
            Z = [Z;gz];
        end
        
        % store
        X = [X;gx];
        T = [T;gt];
        Y = [Y;ysub];
        Y2 = [Y2;y2];
end
%%% concatenate
if isempty(interactionfactor)
    for i = unique(X')
        Y2(X==i) = Y2(X==i) + mean(Y(X==i)); 
    end
else
    [G] = findgroups(X,Z);
    for i = unique(G')
        Y2(G==i) = Y2(G==i) + mean(Y(G==i)); 
    end
end

SUB = repmat([1:nsub],ntrt,1);
SUB = reshape(SUB,nsub*ntrt,1);



% display
clear g;
alpha=0.7;   
if isempty(interactionfactor)
    if sessionSplit
        g = gramm('x',S,'y',Y2,'color',T);
    elseif ~isempty(order)
        g = gramm('x',O,'y',Y2,'color',T);
    else
        g = gramm('x',X,'y',Y2,'color',T);
    end
    g.set_color_options('map',vertcat(col{:}),'lightness',100);
    g.set_order_options('x',treatmentList,'color',treatmentList);
    for i = 1:numel(plotType)
        switch plotType{i}
            case 'dot'
                g.geom_point();
            case 'jitter'
                g.geom_jitter('width',0.05,'dodge',0.3);
            case 'boxplot'
                    g.stat_boxplot('width',0.9);
            case 'errorbar'
                g.stat_summary('type','sem','geom',{'point','line','errorbar'},'width',0.01);
        end
    end
    if sessionSplit==1
        xname = 'session number';
    elseif ~isempty(order)
        xname = 'task order';
    else
        xname = '';
    end
    g.set_names('x',legend{1},'y',legend{2},'color','treatment');
    if sessionSplit || ~isempty(order); g.axe_property('XLim',[0 4]); end
    g.axe_property('YLim',[nanmin(Y2) nanmax(Y2)] + [-0.2 0.2]*nanmean(Y2));
else
    j=0;
    for iz = unique(Z)'
        j=j+1;
        if sessionSplit
            g(1,j) = gramm('x',S,'y',Y2,'color',T);
        elseif ~isempty(order)
            g(1,j) = gramm('x',O,'y',Y2,'color',T);
        else
            g(1,j) = gramm('x',X,'y',Y2,'color',T,'subset',(Z==iz));
        end
        g(1,j).set_color_options('map',vertcat(col{:}),'lightness',100);
        g(1,j).set_order_options('x',treatmentList,'color',treatmentList);
        if isnumeric(iz)
            leg = num2str(iz);
        else
            leg = char(iz);
        end
        g(1,j).set_title([ legend{3} '=' leg  ]);
        for i = 1:numel(plotType)
            switch plotType{i}
                case 'dot'
                    g(1,j).geom_point();
                case 'jitter'
                    g(1,j).geom_jitter('width',0.05,'dodge',0.3);
                case 'boxplot'
                    g(1,j).stat_boxplot('width',0.9);
                case 'errorbar'
                    g(1,j).stat_summary('type','sem','geom',{'point','line','errorbar'},'width',0.01);
            end
        end
        if sessionSplit==1
            xname = 'session number';
        elseif ~isempty(order)
            xname = 'task order';
        else
            xname = '';
        end
        g(1,j).set_names('x',legend{1},'y',legend{2},'color','treatment');
        if sessionSplit || ~isempty(order); g(1,j).axe_property('XLim',[0 4]); end
        g(1,j).axe_property('YLim',[nanmin(Y2) nanmax(Y2)] + [-0.2 0.2]*nanmean(Y2));
    end
end
g.draw;
if isempty(interactionfactor)
    axes(g.facet_axes_handles);
end
set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                    'LineWidth',1.5,'FaceAlpha',alpha);

end