%% analyze_batmotiv_opioid_2level

% reset
clc;
clear all;
close all;

%% Prepare data analysis
%-------------------------------

% load dataset
%%% design table
datadir = 'B:\nicolas.borderies\projets\batmotiv\données\OPIOID';
cd(datadir);
load('subject_table_opioid.mat');
design = tab;  
% processed dataset
datadir = 'B:\nicolas.borderies\projets\batmotiv\resultats\OPIOID';
codedir = 'B:\nicolas.borderies\projets\batmotiv\code.perso';
cd(datadir);
analysisname = 'motiscan-opioid-1level_22_8_2017.mat';
load(analysisname);
cd(codedir);

        
%% Data analysis
%-------------------------------

% define analysis parameters
%%% reformat
data = groupdata;
result = groupresult;
%%% lists
dimensionList = {'writtenReward','visualReward','writtenPunishment','writtenEffort',...
                            'monetaryReward','monetaryPunishment','gripEffort',...
                            'RewardEffort','PunishmentEffort','RewardPunishment',...
                            'Gain','Loss','GainLoss','GainEmotion'};
taskList =  {'rating','choice','weight','discount',...
                    'grip','gripIAPS','gripAccu','mental','learning'};
treatmentList = {'naloxone','placebo','morphine'};
orderGrip = design.order_taskGrip(design.selection==1);
orderRating = design.order_rating(design.selection==1);
orderWeight = design.order_weight(design.selection==1);
weight = design.weight(design.selection==1);

ntrt = numel(treatmentList);
nsub = numel(data);
%%% display
col = { [0 0 1]*0.75 , [1 1 1]*0.5 , [1 0 0]*0.75 };

%% 1) Grip task
%------------------------

taskname = 'grip';

%% Data completion
for isub = 1:nsub % subject loop
        % select
        tab = data{isub}.battery.table; 
        selection =  (tab.task==taskname) ;
        tab = tab(selection,:);
        % variables
        incentive = tab.incentiveLevel ; 
        valence = tab.incentiveSign ; 
        nt = tab.trialNumber ; 
        qnt = quantileranks(nt,2) ; 
        fpeak = tab.normalizedForcePeak ;
        fsum = tab.normalizedForceSum ;
        rt = tab.rt;
        vpeak = tab.normalizedYankPeak;
        trt = tab.treatment;
        trt = removecats(trt,'0');
        session = tab.sessionNumber;
        
        

%         trt = reordercats(trt,treatmentList);
        % stats
            % orthogonalization of force peak & factors
            [g] = findgroups(incentive,qnt);
            [~,~,ycond] = varpart(fpeak,g);
            data{isub}.(taskname).table.cond_normalizedForcePeak = ycond;
            [g] = findgroups(trt);
            mean_cond_normFPeak = splitapply(@nanmean,ycond,g);
            result{isub}.(taskname).inferential.mean_cond_normFPeak = mean_cond_normFPeak;
            
            % glm of force peak
            intercept = nan(3,1);
            for isess=1:3
                predictor = [incentive,incentive.*valence,nt];
                beta = glmfit(predictor(session==isess,:),fpeak(session==isess),'normal');
                intercept(isess) = beta(1);
            end
            result{isub}.(taskname).inferential.intercept_norm_fpeak = intercept;       
            
            % discretization
            discrete_fpeak = quantileranks(fpeak,6);
            discrete_fpeak(discrete_fpeak==0)=NaN;
            data{isub}.(taskname).table.discrete_fpeak = discrete_fpeak;       
            
            % force peak - velocity 
            intercept = nan(3,1);
            for isess=1:3
                beta = glmfit(fpeak(session==isess,:),vpeak(session==isess),'normal');
                intercept(isess) = beta(2);
            end
            result{isub}.(taskname).inferential.beta_force_velocity = intercept;  

            % force peak - rt
            intercept = nan(3,1);
            for isess=1:3
                beta = glmfit((fpeak(session==isess,:)),rt(session==isess),'normal');
                intercept(isess) = beta(2);
            end
            result{isub}.(taskname).inferential.beta_rt_force = intercept;  
            
            data{isub}.(taskname).table.fpeakByTime = fpeak./rt;
            data{isub}.(taskname).table.fsumByTime = fsum./rt;
            data{isub}.(taskname).table.fmean = cellfun(@nanmean,data{isub}.grip.table.normalizedForceValue(:));
            data{isub}.(taskname).table.fsum = cellfun(@nanmean,data{isub}.grip.table.normalizedForceValue(:));

            
end

%% Univariate effect of treatments:
% yname = 'normalizedForcePeak';
yname = 'normalizedForceSum';
% yname = 'cond_normalizedForcePeak';
% yname = 'trialGain';
% yname = 'time2forcePeak';
% yname = 'rt';
% yname = 'normalizedYankPeak';
% yname = 'fpeakByTime';
% yname = 'fsumByTime';
% yname = 'fmean';


% ylegend = 'force peak (%fmax)';
% ylegend = 'velocity peak (%fmax/sec)';
% ylegend = 'maximal voluntary force peak (%fmax)';
% ylegend = 'preparation time (sec)';
ylegend = ' cumulative effort (%fmax)';


fname = @nanmean;
% fname = @nanmax;
% ftransform = @identity;
% ftransform = @ (x) normalize(x,'zscore');
ftransform = @ (x) normalize(x,'max');

sessionEffect=0;
%   - statistic:
statistic = {'opioid_ttest'};
% statistic = {'naloxone_contrast','morphine_contrast'};
% statistic = {'naloxone_ttest','morphine_ttest'};
%   - plot: 'bar','dot','jitter','boxplot','violin','line'
plot = {'jitter','boxplot'};

[p,score] = barcomp_motiscan_opioid(data,yname,ylegend,taskname,...
                               'fname',fname,'ftransform',ftransform,...
                               'sessionSplit',sessionEffect,...
                               'statisticalTest',statistic,...
                               'plotType',plot); % without session interaction

%% Parametric effect of treatments:

%   - predicted variables
% yname = {'offset_normalizedForcePeak','incentive_normalizedForcePeak','incentive_gain_normalizedForcePeak','trialNumber_normalizedForcePeak'};
% yname = {'kR','kP','kE','kF','tau'};
% yname = {'fmax'};
% yname = {'offset_normalizedForcePeak'};
% yname = {'intercept_norm_fpeak'};
% yname = {'beta_force_velocity'};
yname = {'beta_rt_force'};

%   - transform
% ftransform = @log;
ftransform = @identity;
% ftransform = @ (x) normalize(x,'variance');
%   - legend
ylegend = '';
% ylegend = 'parameters';
%   - order interaction
order = [];
%   - statistic:
% statistic = {'opioid_ttest'};
statistic = {'naloxone_contrast','morphine_contrast'};
% statistic = {'naloxone_ttest','morphine_ttest'};
%   - plot: 'bar','dot','jitter','boxplot','violin','line'
plot = {'boxplot'};

[p,score,stat] = barcomp_stat_motiscan_opioid(result,yname,ylegend,taskname,...
                               'fname',fname,...
                               'ftransform',ftransform,...
                               'sessionSplit',sessionEffect,...
                               'order',order,...
                               'statisticalTest',statistic,...
                               'plotType',plot);
                           
%% Bivariate effects of treatment 


%   - predicted variable
% yname = 'normalizedForcePeak';
yname = 'normalizedForceSum';
% yname = 'rt';
% yname = 'normalizedYankPeak';
%   - factor
xname = 'incentiveLevel';
% xname = 'incentiveSign';
% xname = 'trialNumber';
% xname = 'discrete_fpeak';
interaction = [];


%   - transform
ftransform = @identity;
ftransform = @ (x) normalize(x,'max');
%   - statistical function
fname = @nanmean;
% fname = @nanmax;
%   - session interaction
sessionEffect=0;
%   - legend
legend = {'incentive','cumulative effort (%fmax)'};
%   - order interaction
order = [];
% order = orderGrip;
%   - plot: 'dot','jitter','line'
% plot = {'dot','errorbar'};
% plot = {'jitter','errorbar'};
plot = {'errorbar'};

figure;
bivarcomp_motiscan_opioid(data,xname,yname,legend,taskname,...
                               'fname',fname,...
                               'ftransform',ftransform,...
                               'sessionSplit',sessionEffect,...
                               'order',order,...
                               'statisticalTest',statistic,...
                               'plotType',plot,...
                               'interactionfactor',interaction);
                           
%% Treatments-Cofactors interactions
yname = 'normalizedForcePeak';
factor = weight;
legend = {'weight','force peak (%fmax)'};

ftransform = @identity;
fname = @nanmean;
% fname = @nanmax;
%   - plot: 'dot','jitter','line'
plot = {'fit','jitter'};

figure;
cofactor_comp_opioid(data,factor,yname,legend,taskname,...
                               'fname',fname,...
                               'ftransform',ftransform,...
                               'plotType',plot);
                           
                           
%% Parametric effect of treatments (GLM)
%   - predicted variable
% yname = 'normalizedForcePeak';
yname = 'normalizedForceSum';
% yname = 'rt';
%   - model equation
formula = [ yname ' ~ 1 + incentive + valence:incentive + ntrial',...
           '+ naloxone:(1 + incentive + valence:incentive + ntrial)',...
           '+ morphine:(1 + incentive + valence:incentive + ntrial)',...
           '+ session'];
% formula = [ yname ' ~ 1 + incentive + fpeak + valence:incentive',...
%            '+ naloxone:(1 + incentive + fpeak + valence:incentive)',...
%            '+ morphine:(1 + incentive + fpeak + valence:incentive )',...
%            '+ session'];
       
varnames = {'incentive','valence','ntrial','fpeak',...
            'naloxone','placebo','morphine','opioid','session',yname};
        
ftransform = @identity;
% ftransform = @log;
%   - statistic:
%   'drug_ttest','opioid_ttest'
statistic = {'opioid_ttest'};
% statistic = {'drug_ttest'};
%   - plot: 'bar','dot','jitter','boxplot','violin','line'
plot = {'boxplot'};
% plot = {'bar'};


figure;
[p,coefNames,stat] = glm_grip_motiscan_opioid(data,yname,varnames,formula,...
                                                         'plotType',plot,...
                                                         'statisticalTest',statistic,...
                                                         'ftransform',ftransform);
                           
%%
% raw effect

    nbin = [ 3 , nsub ];

    % variables
        Y = nan(nbin);
        Y2 = nan(nbin);


    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='grip') ;
            tab = tab(selection,:);
            
        % variables
%             force = tab.forcePeak;
            force = tab.normalizedForcePeak;
%             force = tab.normalizedForceSum;
%             force = tab.rt;
%             force = tab.yankPeak;

             trt = tab.treatment;
             trt = removecats(trt,'0');
             session = tab.sessionNumber;
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            trt = reordercats(trt,treatmentList(subtrt));
            
            ysub = tools.tapply(force,{trt},@nanmean);
%             ysub = tools.tapply(force,{trt},@nanmax);

            Y(subtrt,isub) =  ysub;
            Y2(subtrt,isub) =  Y(subtrt,isub)- Y(2,isub);

    end
    
    % display
        fig = figure; set(fig,'Name','grip_opioid');
            
            dsub=2;
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);
            
            hold on
            for it=1:3
                x = it;
                [ h(it)  ] = barplot( x ,y(it),z(it), col{it} );
                h(it).BarWidth = 0.5;
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            ax = gca; 
            ax.XTick = [];
%             ylabel('force peak (%fmax) '); 
            ylabel('maximal force peak (N) '); 

            ax.XLim = [0 x(end)+1];
            ax.YLim = [min(y) max(y)] + [-0.2 0.2]*mean(y);

        % format
        setFigProper('FontSize',20,'LineWidth',2);


%% second-level parameters stats
% GLM  parameters

  nbin = [ 3 , 3, nsub ];

    % variables
        Y = nan(nbin);
        Y2 = nan(nbin);

    for isub = 1:nsub
%         try
        % select
            tab = result{isub}.grip.inferential  ; 
            
        % variables
%             kI = nanmean([tab.kR tab.kP],2);
            kI = tab.incentive_forcePeak ;
            kIG = tab.incentive_gain_forcePeak ;
            kI = (kI + kIG)/2 ;

            kE = tab.offset_forcePeak;
            kE = 1 - tab.mean_normalizedForcePeak;
            kF = tab.trialNumber_forcePeak ;

            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            session = [1 2 3]';

            
        % stats
            [~,subtrt] = ismember(trt,treatmentList);
            ntrt = numel(subtrt);
            ysub = [kE kI kF ];
            Y(subtrt,:,isub) =  ysub;
            Y2(subtrt,:,isub) =  Y(subtrt,:,isub)- repmat(Y(2,:,isub),ntrt,1);
%         end

    end
    
    % display
        fig = figure; set(fig,'Name','grip_param');
            
        % stat
        dsub=3;
        y = nanmean(Y,dsub);
        z = sem(Y2,dsub);     
        % barplot
        hold on; clear h ; 
        for it=1:3
            x = [1 2 3];
            xx = x + (it-2)*0.2;
            [ h(it)  ] = barplot( xx ,y(it,:),z(it,:), col{it} );
            %                 [ h(it) ] = myboxplot( x , reshape(Y(it,:),nsub,1) , col{it} );
            h(it).BarWidth = 0.2;
            
            [hyp,p] = ttest(reshape(Y2(it,:,:),nbin(2),nsub)');
            [s] = sigstar( num2cell(xx),p);
        end
        legend([h(1) h(2) h(3)],treatmentList);

        % legending
        paramName = {'kE','kI','kF'};
        ax = gca; 
        ax.XTick = x;
        ax.XTickLabel = paramName(x);
        ylabel(' parametric effects  '); 
        ax.XLim = [0 max(x)+1];

%% second-level parameters stats
% nonlinear model parameters
  nparam = 5;
  nbin = [ 3 , nparam , nsub ];

    % variables
        Y = nan(nbin);

    for isub = 1:nsub
        try
        % select
            tab = result{isub}.grip.inferential  ; 
            
        % variables

            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            session = [1 2 3]';

            
        % stats
            [~,subtrt] = ismember(trt,treatmentList);
            ysub = [tab.kR  tab.kP  tab.kE tab.kF tab.tau ];

            Y(subtrt,:,isub) =  ysub;
        end

    end
    
    % display
        fig = figure; set(fig,'Name','pref_param');
            
            dsub=3;
            y = nanmedian(Y,dsub);
            z = sem(Y,dsub);
            ylist = {'kR','kP','kE','kF','tau'};
            
            for idim=1:nparam
                subplot(1,nparam,idim);
                hold on;
                for it=1:3
                    x = it;
                    [ h(it)  ] = barplot( x ,y(it,idim),z(it,idim), col{it} );
    %                 [ h(it) ] = myboxplot( x , reshape(Y(it,:),nsub,1) , col{it} );
                    h(it).BarWidth = 0.5;
                end
                
                legend([h(1) h(2) h(3)],treatmentList);            
                % legending
                ax = gca; 
                ax.XTick = [];
                title(ylist{idim}); 
            end
        
                
%% incentive effect

    nbin = [ 3 , 6 , nsub ];

    % variables
        X = nan(nbin);
        Y = nan(nbin);
        Y2 = nan(nbin);

    for isub = 1:nsub
%         try
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='grip') ;
            tab = tab(selection,:);
            
        % variables
            force = tab.normalizedForcePeak;
%             force = tab.normalizedForceSum;
%             force = tab.rt;
%             prediction = tab.predicted_forcePeak./tab.maxObservedForcePeak;


            incentive = tab.incentiveValue;
            valence = tab.incentiveSign;
            
            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            X(subtrt,:,isub) =  tools.tapply(incentive,{trt,incentive},@nanmean);
            ysub =  tools.tapply(force,{trt,incentive},@nanmean);
            ysub =  tools.tapply(force,{trt,incentive,valence},@nanmean);
            Y(subtrt,:,isub) = ysub(:,:,1);
            
%             Y2(subtrt,:,isub) =  tools.tapply(prediction,{trt,incentive},@nanmean);

%         end

    end
    
    % display            
            dsub=3;
            x = nanmean(X,dsub);
            y = nanmean(Y,dsub);
%             y2 = nanmean(Y2,dsub);
            z = sem(Y,dsub);
            
            % group
                fig = figure; set(fig,'Name','grip_opioid');
                hold on
                for it=1:3
                    x = [1:6];
                    [~,~,h(it)] = errorscat(x,y(it,:),z(it,:),col{it});
%                     h(it).LineStyle='none';
%                     h(it) = plot(x,y2(it,:),'Color',col{it});
                end
                legend([h(1) h(2) h(3)],treatmentList);
            
            % subject-wise
%                 for isub = 1:nsub
%                     if exist('fig');close(fig); clear fig ; end
%                     fig = figure; hold on; 
%                     x = X(:,:,isub);
%                     y = Y(:,:,isub);
%                     y2 = Y2(:,:,isub);
% 
%                     for it=1:3
%                         x = [1:6];
%                         [~,~,h(it)] = errorscat(x,y(it,:),y(it,:).*0,col{it});
%                         h(it).LineStyle='none';
%                         h(it) = plot(x,y2(it,:),'Color',col{it});
%                     end
%                 legend([h(1) h(2) h(3)],treatmentList);
%                 pause;
%                 end
            
        % legending
            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = cellfun(@num2str,num2cell(unique(incentive)),'UniformOutput',0);
            xlabel('incentive(€)'); 
            ylabel('force peak (%fmax) '); 
            ax.XLim = [0 x(end)+1];
%             ax.YLim = [0 1];

%% fatigue effect
    nbin = [ 3 , 10 , nsub ];

    % variables
        X = nan(nbin);
        Y = nan(nbin);

    for isub = 1:nsub
        try
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='grip') ;
            tab = tab(selection,:);
            
        % variables
            force = tab.normalizedForcePeak;
%             force = tab.normalizedForceSum;
%             force = tab.rt;

            incentive = tab.incentiveValue;    
            nt = tab.trialNumber;
            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            X(subtrt,:,isub) =  tools.tapply(nt,{trt,nt},@nanmean,{'discrete','continuous'},nbin(1:2));
            Y(subtrt,:,isub) =  tools.tapply(force,{trt,nt},@nanmean,{'discrete','continuous'},nbin(1:2));
        end
    end
    
    % display
        fig = figure; set(fig,'Name','grip_opioid');
            
            dsub=3;
            x = nanmean(X,dsub);
            y = nanmean(Y,dsub);
            z = sem(Y,dsub);
            

            hold on
            for it=1:3
                [~,~,h(it)] = errorscat(x(it,:),y(it,:),z(it,:),col{it});
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            ax = gca; 
%             ax.XTick = x(it,:)-0.5;
            ax.XTick = [10:10:120];

            xlabel('trial number'); 
            ylabel('force peak (%fmax) '); 
            ax.XLim = [0 x(end)+10];
%             ax.YLim = [0 1];

        
%% 2) Sustained Effort task
%------------------------------------------

taskname = 'gripAccu';

%% Data completion
for isub = 1:nsub % subject loop
        % select
        tab = data{isub}.battery.table; 
        selection =  (tab.task==taskname) ;
        tab = tab(selection,:);
        % variables
        y = tab.forceDuration ./(tab.trialDuration);
%         data{isub}.battery.table.effortRatio = [];
        data{isub}.(taskname).table{:,'effortRatio'} = y;
        y = tab.trialDuration -  tab.responseTime;
        data{isub}.(taskname).table{:,'squeezeDuration'} = y;
        meanTD = nanmean(tab.trialDuration);
        maxTD = nanmax(tab.trialDuration);
        data{isub}.(taskname).table{:,'prepareDuration_meanTD'} = tab.preparatoryDuration./meanTD;
        data{isub}.(taskname).table{:,'prepareDuration_maxTD'} = tab.preparatoryDuration./maxTD;
end

%% Univariate effect of treatments:

%   - predicted variable
% yname = 'effortDuration';
% yname = 'restDuration';
yname = 'effortRatio';
% yname = 'forceDuration';
% yname = 'preparatoryDuration';
% yname = 'prepareDuration_meanTD';
% yname = 'prepareDuration_maxTD';
% yname = 'norm_forceAmplitude';
% yname = 'correct_force';
% yname = 'norm_forceSum';
% yname = 'responseTime';
% yname = 'squeezeDuration';
% yname = 'gain';
% yname = 'trialDuration';
%   - transform
% ftransform = @ (x) normalize(x,'zscore');
ftransform = @identity;
%   - statistical function
fname = @nanmean;
% fname = @nanmax;
%   - session interaction
sessionEffect=0;
%   - legend
% ylegend = 'force duration (sec)';
% ylegend = 'preparation duration (sec)';
ylegend = 'effort rate (%)';

%   - order interaction
order = [];
% order = orderGrip;
%   - statistic:
%   'morphine_contrast','naloxone_contrast','opioid_contrast',...
%   'morphine_ttest','naloxone_ttest','opioid_ttest'
statistic = {'opioid_contrast'};
% statistic = {'naloxone_contrast','morphine_contrast'};
% statistic = {'naloxone_ttest','morphine_ttest'};
%   - plot: 'bar','dot','jitter','boxplot','violin','line'
plotType = {'jitter','boxplot'};
% plotType = {'boxplot','line','dot'};


figure;
[p,score,g] = barcomp_motiscan_opioid(data,yname,ylegend,taskname,...
                               'fname',fname,...
                               'ftransform',ftransform,...
                               'sessionSplit',sessionEffect,...
                               'order',order,...
                               'statisticalTest',statistic,...
                               'plotType',plotType);
                           
%% Bivariate effects of treatment 


%   - predicted variable
% yname = 'effortDuration';
yname = 'effortRatio';
% yname = 'forceDuration';
% yname = 'preparatoryDuration';
% yname = 'norm_forceAmplitude';
% yname = 'correct_force';
% yname = 'norm_forceSum';
% yname = 'responseTime';
% yname = 'squeezeDuration';
% yname = 'gain';
% yname = 'trialDuration';
%   - factor
% xname = 'incentiveLevel';
xname = 'costLevel';
interaction = 'incentiveLevel';
%   - transform
ftransform = @identity;
%   - statistical function
fname = @nanmean;
% fname = @nanmax;
%   - session interaction
sessionEffect=0;
%   - legend
legend = {'costLevel','force duration (% trial duration)','incentive'};
%   - order interaction
order = [];
% order = orderGrip;
%   - plot: 'dot','jitter','line'
% plot = {'dot','errorbar'};
% plot = {'jitter','errorbar'};
plot = {'errorbar'};

%   - interaction


figure;
bivarcomp_motiscan_opioid(data,xname,yname,legend,taskname,...
                               'fname',fname,...
                               'ftransform',ftransform,...
                               'sessionSplit',sessionEffect,...
                               'order',order,...
                               'statisticalTest',statistic,...
                               'plotType',plot,...
                               'interactionfactor',interaction);

%% Treatments-Cofactors interactions
yname = 'effortRatio';
% yname = 'forceDuration';
factor = weight;
legend = {'weight','force duration (%trial duration)'};

ftransform = @identity;
fname = @nanmean;
% fname = @nanmax;
%   - plot: 'dot','jitter','line'
plot = {'fit','jitter'};

figure;
cofactor_comp_opioid(data,factor,yname,legend,taskname,...
                               'fname',fname,...
                               'ftransform',ftransform,...
                               'plotType',plot);
                           
%% Parametric effect of treatments (GLM)
%   - predicted variable
% yname = 'forceDuration';
% yname = 'preparatoryDuration';
yname = 'effortRatio';
%   - model equation
formula = [ yname ' ~ 1 + incentive + difficulty + instruction:difficulty + nblock + ntrial_block + hand',...
           '+ naloxone:(1 + incentive + difficulty + instruction:difficulty + nblock + ntrial_block + hand)',...
           '+ morphine:(1 + incentive + difficulty + instruction:difficulty + nblock + ntrial_block + hand)',...
           '+ session'];
% formula = [ yname ' ~ 1 + incentive + difficulty + instruction:difficulty + nblock + ntrial_block + hand',...
%            '+ opioid:(1 + incentive + difficulty + instruction:difficulty + nblock + ntrial_block + hand)',...
%            '+ session'];
varnames = {'incentive','difficulty','instruction','effort_t',...
            'nblock','ntrial_block','hand',...
            'naloxone','placebo','morphine','opioid','session',yname};
%   - statistic:
%   'drug_ttest','opioid_ttest'
% statistic = {'opioid_ttest'};
statistic = {'drug_ttest'};
%   - plot: 'bar','dot','jitter','boxplot','violin','line'
plot = {'boxplot'};
% plot = {'bar'};

figure;
[p,coefNames,stat] = glm_gripAccu_motiscan_opioid(data,yname,varnames,formula,...
                                                         'plotType',plot,...
                                                         'statisticalTest',statistic);
                    
                    
                    
%% cost effect
    nbin = [ 3 , 2 , nsub ];

    % variables
        Y = nan(nbin);
        Y2 = nan(nbin);

        
    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            tab = data{isub}.('gripAccu').table; 
            selection =  (tab.task=='gripAccu') ;
            tab = tab(selection,:);
            
        % variables
            effort = tab.effortDuration;
%             effort = tab.predicted_effortDuration;

            rest = tab.restDuration;

            cost = tab.costValue;
            incentive = tab.incentiveValue;

            trt = tab.treatment;
            trt = removecats(trt,'0');
            session = tab.sessionNumber;
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            trt = reordercats(trt,treatmentList(subtrt));
            ntrt = numel(subtrt);
            
%             ysub = tools.tapply(effort,{trt,cost},@nanmean);
            ysub = tools.tapply(rest,{trt,cost},@nanmean);

            
            Y(subtrt,:,isub) =  ysub;
            Y2(subtrt,:,isub) =  Y(subtrt,:,isub)- repmat(Y(2,:,isub),ntrt,1);


    end
    
    % display
        fig = figure; set(fig,'Name','gripAccu_opioid_cost');
            
            dsub=3;
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);

            hold on;
            for it=1:3
                x = [1:2];
                [~,~,h(it)] = errorscat(x,y(it,:),z(it,:),col{it});
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = {'75%','85%'};
            xlabel('difficulty level (%fmax)  '); 
            ylabel('force duration (sec) '); 
%             ylabel('rest duration (sec) '); 

            ax.XLim = [0 x(end)+1];
%             ax.YLim = [0 1];


        % format
        setFigProper('FontSize',20,'LineWidth',2);

%% instruction effect
    nbin = [ 3 , 2 , 2 , nsub ];

    % variables
        Y = nan(nbin);

    for isub = 1:nsub
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='gripAccu') ;
            tab = tab(selection,:);
            
        % variables

            effort = tab.effortDuration;
%             effort = tab.predicted_effortDuration;
%             rest = tab.restDuration;
%             rest = tab.predicted_restDuration;

            cost = tab.costValue;
            incentive = tab.incentiveValue;
            instruction = tab.explicitCost;
            
            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            session = tab.sessionNumber;
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
%             [~,subsess] = ismember(unique(session),[1 2 3]);
            Y(subtrt,:,:,isub) =  tools.tapply(effort,{trt,cost,instruction},@nanmean);


    end
    
    % display
        fig = figure; set(fig,'Name','gripAccu_opioid_instruction');
            
            dsub=4;
            y = nanmedian(Y,dsub);
            z = sem(Y,dsub);

            hold on;
            for it=1:3
                for ins=1:2
                    x = [1:2];
                    [~,~,h(it)] = errorscat(x,y(it,:,ins,1),z(it,:,ins,1),col{it});
                    if ins==1; h(it).LineStyle = '--'; end
                end
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = {'75%','85%'};
            xlabel('force level (%fmax)  '); 
            ylabel('force duration (zscore|subject) '); 
            ax.XLim = [0 x(end)+1];
%             ax.YLim = [0 1];

%% incentive effect
     nbin = [ 3 , 2 , nsub ];

    % variables
        Y = nan(nbin);
        Y2 = nan(nbin);

        
    for isub = 1:nsub
                
        % select
%             tab = data{isub}.battery.table; 
            tab = data{isub}.('gripAccu').table; 
            selection =  (tab.task=='gripAccu') ;
            tab = tab(selection,:);
            
        % variables

%             effort = tab.effortDuration;
            effort = tab.effortRatio;

%             effort = tab.predicted_effortDuration;
            rest = tab.restDuration;

            cost = tab.costValue;
            incentive = tab.incentiveValue;

            trt = tab.treatment;
            trt = removecats(trt,'0');
            session = tab.sessionNumber;
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            trt = reordercats(trt,treatmentList(subtrt));
            ntrt = numel(subtrt);
            
            ysub = tools.tapply(effort,{trt,incentive},@nanmean);
%             ysub = tools.tapply(rest,{trt,incentive},@nanmean);

            
            Y(subtrt,:,isub) =  ysub;
            Y2(subtrt,:,isub) =  Y(subtrt,:,isub)- repmat(Y(2,:,isub),ntrt,1);


    end
    
    % display
        fig = figure; set(fig,'Name','gripAccu_opioid_incentive');
            
            dsub=3;
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);

            hold on;
            for it=1:3
                x = [1:2];
                [~,~,h(it)] = errorscat(x,y(it,:),z(it,:),col{it});
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = {'0.10','0.20'};
            xlabel('incentive level (€)  '); 
            ylabel('force duration (sec) '); 
            ylabel('effort rate (%fmax) '); 
%             ylabel('rest duration (sec) '); 

            ax.XLim = [0 x(end)+1];
%             ax.YLim = [0 1];


        % format
        setFigProper('FontSize',20,'LineWidth',2);


%% previous trial effect
    nbin = [ 3 , 6 , nsub ];

    % variables
        X = nan(nbin);
        Y = nan(nbin);
        Y2 = nan(nbin);

        
    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='gripAccu') ;
            tab = tab(selection,:);
            
        % variables
            effort = tab.effortDuration;
            rest = tab.restDuration;

            cost = tab.costValue;
            incentive = tab.incentiveValue;
            nt = tab.trialNumber;
            
            effort_t = nan(size(nt));
            effort_t(2:end) = effort(1:end-1);
            effort_t(nt<=1) = NaN;
            cost_t = nan(size(nt));
            cost_t(2:end) = cost(1:end-1);
            cost_t(nt<=1) = NaN;
            
            cumCost = cost_t.*effort_t;
            cumCost = effort_t;
            cumCost2 = quantileranks(cumCost,nbin(2));
            cumCost2(cumCost2==0)=1;
            
            trt = tab.treatment;
            trt = removecats(trt,'0');
            session = tab.sessionNumber;
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            trt = reordercats(trt,treatmentList(subtrt));
            ntrt = numel(subtrt);
            
            xsub = tools.tapply(cumCost,{trt,cumCost2},@nanmean);
            ysub = tools.tapply(rest,{trt,cumCost2},@nanmean);
            
            X(subtrt,:,isub) =  xsub;
            Y(subtrt,:,isub) =  ysub;
            Y2(subtrt,:,isub) =  Y(subtrt,:,isub)- repmat(Y(2,:,isub),ntrt,1);

    end
    
    % display
        fig = figure; set(fig,'Name','gripAccu_opioid_incentive');
            
            dsub=3;
            x = nanmean(X,dsub);
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);

            hold on;
            for it=1:3
%                 x = [1:nbin(2)];
                [~,~,h(it)] = errorscat(x(it,:),y(it,:),z(it,:),col{it});
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            ax = gca; 
%             ax.XTick = x;
%             ax.XTickLabel = {'0.10','0.20'};
            xlabel('force duration_{t-1}  (sec) '); 
            ylabel('rest duration (sec) '); 
%             ax.XLim = [0 max(max(x))+1];
%             ax.YLim = [0 1];


    % format
    setFigProper('FontSize',20,'LineWidth',2);

%% fatigue effect
    nbin = [ 3 , 32 , nsub ];

    % variables
        X = nan(nbin);
        Y = nan(nbin);
        Y2 = nan(nbin);

        
    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='gripAccu') ;
            tab = tab(selection,:);
            
        % variables
            effort = tab.effortDuration;
            rest = tab.restDuration;

            cost = tab.costValue;
            incentive = tab.incentiveValue;
            nt = tab.trialNumber;
           
            trt = tab.treatment;
            trt = removecats(trt,'0');
            session = tab.sessionNumber;
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            trt = reordercats(trt,treatmentList(subtrt));
            ntrt = numel(subtrt);
            
            xsub = tools.tapply(nt,{trt,nt},@nanmean);
%             ysub = tools.tapply(effort,{trt,nt},@nanmean);
            ysub = tools.tapply(rest,{trt,nt},@nanmean);

            X(subtrt,:,isub) =  xsub;
            Y(subtrt,:,isub) =  ysub;
            Y2(subtrt,:,isub) =  Y(subtrt,:,isub)- repmat(Y(2,:,isub),ntrt,1);

    end
    
    % display
        fig = figure; set(fig,'Name','gripAccu_opioid_fatigue');
            
            dsub=3;
            x = nanmean(X,dsub);
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);

            hold on;
            for it=1:3
%                 x = [1:nbin(2)];
                [~,~,h(it)] = errorscat(x(it,:),y(it,:),z(it,:),col{it});
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            ax = gca; 
%             ax.XTick = x;
%             ax.XTickLabel = {'0.10','0.20'};
            xlabel('n(trial) '); 
%             ylabel('effort duration (sec) '); 
            ylabel('rest duration (sec) '); 

%             ax.XLim = [0 max(max(x))+1];
%             ax.YLim = [0 1];


    % format
    setFigProper('FontSize',20,'LineWidth',2);
    



%% second-level stats
% regressions
%     nbin = [ 3 , 8 , nsub ];
    nbin = [ 3 , 7 , nsub ];
    nbin2 = [ 3 , nsub ];

 % variables
    Y = nan(nbin);
    Y2 = nan(nbin);
    Y3 = nan(nbin2);

    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='gripAccu') ;
            tab = tab(selection,:);
            
        % variables
            effort = tab.effortDuration;
            rest = tab.restDuration;
%             effort = errorscore(tab.effortDuration);
%             rest = errorscore(tab.restDuration);

            cost = tab.costValue;
            incentive = tab.incentiveValue;
            instruction = tab.explicitCost;
            
            trial_task = tab.trialNumber;
            trial_hand = mod(trial_task-1,8)+1;
            hand = tab.handSide;
            
            effort_t = nan(size(trial_task));
            effort_t(2:end) = effort(1:end-1);
%             effort_t(2:end) = cost(1:end-1);
            effort_t(trial_task<=1) = NaN;
            
            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
           
            session = tab.sessionNumber;
            
        % first level stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            for it = subtrt'
                
                predictor = ([incentive,cost,instruction,effort_t,trial_task,trial_hand,hand]);
                predictor = nanzscore(predictor(trt==treatmentList(it),:));
                indcat =  [3];
                predictor(:,indcat) = predictor(:,indcat).*(predictor(:,indcat)>0) ;

                y = effort(trt==treatmentList(it));
%                 y = 1./effort(trt==treatmentList(it));
                formula = 'effort ~ 1 + incentive + cost + instruction:cost + trial_task + trial_hand + hand';
                varnames = {'incentive','cost','instruction','effort_t','trial_task','trial_hand','hand','effort'};

%                 y = rest(trt==treatmentList(it));
%                 formula = 'rest ~ 1 + incentive + effort_t + instruction:cost + trial_task + trial_hand + hand ';
%                 varnames = {'incentive','cost','instruction','effort_t','trial_task','trial_hand','hand','rest'};
%                 
                stat = fitglm(predictor,y,...
                               formula,...
                               'VarNames',varnames);
                coef = stat.Coefficients;
                disp(coef);
               Y(it,:,isub) = coef.Estimate;
               
               Y3(it,isub) = unique(session(trt==treatmentList(it)));

            end
            % average normalization
%             Y(:,:,isub) = errorscore(Y(:,:,isub));
%             denominator = abs(Y(2,:,isub));
%             Y(:,:,isub) = Y(:,:,isub)./repmat(denominator,3,1);
%             denominator = abs(Y(2,1,isub));
%             Y(:,:,isub) = Y(:,:,isub)./denominator;
            comparator = Y(2,:,isub);
            Y2(:,:,isub) =  Y(:,:,isub)- repmat(comparator,3,1);
    end
    
    %% classification
    features = permute(Y,[2 1 3]);
    features = reshape(features,[nbin(2) 24*3])';
    treatment = ordinal(repmat([1;2;3],24,1));
    session = ordinal(reshape(Y3,[1 24*3])');

    [beta,~,stat] = mnrfit(features,treatment,'model','ordinal');
    [post] = mnrval(beta,features,stat,'model','ordinal');
    [~,cat] = max(post,[],2);
    accuracy = mean(treatment==ordinal(cat));
    bca = nan(3,1);
    for it=1:3
        bca(it) = mean(cat(double(treatment)==it)==it);
    end
    accuracy = mean(bca([1 3]));
    chance = 0.33;
    
    %% regression
%     features = permute(Y,[2 1 3]);
    features = permute(Y2,[2 1 3]);
    features = reshape(features,[nbin(2) 24*3])';
    treatment = (repmat([1;2;3],24,1));
    session = (reshape(Y3,[1 24*3])');
    
    y = nan(3,size(features,2));
    z = nan(3,size(features,2));
    p = nan(3,size(features,2));

    for i = 1:size(features,2)
        predictor = [(treatment==1),(treatment==3),session];
        [beta,~,stat] = glmfit(predictor,features(:,i),'normal');
        y(:,i) = beta(2:4);
        z(:,i) = stat.se(2:4);
        p(:,i) = stat.p(2:4);

    end
    
    % display regression
        fig = figure; set(fig,'Name','regression_gripAccu');
             
        predictors = {'naloxone','morphine','session'};
        effectName = {'intercept','incentive','cost','trial_{task}','trial_{hand}','hand','cost_{explicit}'};
%         effectName = {'intercept','incentive','effort_{t-1}','trial_{task}','trial_{hand}','hand','cost_{explicit}'};

        hold on; clear h;
        for i = 1:3
             x = [1:size(features,2)];
             xx = x + (i-2)*0.25 ;
             [ h(i)  ] = barplot( xx ,y(i,:),z(i,:), col{2} );
             h(i).BarWidth = 0.2;
             [s] = sigstar( num2cell(xx),p(i,:));
        end
        h(1).FaceColor = col{1};
        h(2).FaceColor = col{3};

        % legending
            legend([h(1) h(2) h(3)],predictors)
            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = effectName(x);
            ylabel(' coefficients'); 
            ax.XLim = [0 x(end)+1];
    
    %% display glm
        fig = figure; set(fig,'Name','gripAccu_opioid');
            
            dsub=3;
%             y = nanmedian(Y,dsub);
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);
            
            hold on; clear h;
            for it=1:3
%                 x = [1 2 3 4 5 6 7 8];
                x = [1:nbin(2)];

                xx = x + (it-2)*0.2;
                [ h(it)  ] = barplot( xx ,y(it,:),z(it,:), col{it} );
% %                 [ h(it)  ] = barplot( x ,y(it),z(it), col{3} );
                h(it).BarWidth = 0.2;
                
%                 [ h ] = plotSpread(Y','distributionColors',col);
%                 m = findobj('-property','MarkerFaceColor');
%                 set(m,'MarkerSize',14);
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
        effectName = {'intercept','incentive','cost','trial_{task}','trial_{hand}','hand','cost_{explicit}'};
%         effectName = {'intercept','incentive','effort_{t-1}','trial_{task}','trial_{hand}','hand','cost_{explicit}'};


            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = effectName(x);
            ylabel('regression coefficients (sec/var unit) '); 
            title('GLM of effort duration ');  
%             title('GLM of rest duration '); 
            ax.XLim = [0 x(end)+1];
            
%% display classification
        fig = figure; set(fig,'Name','classification_gripAccu');
            
             y = -beta(3:end);
             z = stat.se(3:end);
             p = stat.p(3:end);
             
             hold on; clear h;
             x = [1:nbin(2)];
             xx = x ;
             [ h(it)  ] = barplot( xx ,y,z, col{2} );
             h(it).BarWidth = 0.2;
             [s] = sigstar( num2cell(x),p);
                
        % legending
            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = effectName(x);
            ylabel('regression coefficients (sec/var unit) '); 
            title('treatment classification '); 
            ax.XLim = [0 x(end)+1];
            
    
            
            
    % format
    setFigProper('FontSize',20,'LineWidth',2);

%% gripAccu model
% regressions
    nbin = [ 3 , 10 , nsub ];

 % variables
    Y = nan(nbin);
    Y2 = nan(nbin);

    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='gripAccu') ;
            tab = tab(selection,:);
            
        % variables
            effort = tab.effortDuration;
            rest = tab.restDuration;

            cost = tab.costValue;
            incentive = tab.incentiveValue;
            instruction = tab.explicitCost;
            trial_task = tab.trialNumber;
            trial_hand = mod(trial_task-1,8)+1;
            hand = tab.handSide;
            
            effort_t = nan(size(trial_task));
            effort_t(2:end) = effort(1:end-1);
            effort_t(trial_task<=1) = NaN;
            cost_t = nan(size(trial_task));
            cost_t(2:end) = cost(1:end-1);
            cost_t(trial_task<=1) = NaN;

            
            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
           
            session = tab.sessionNumber;
            
        % first level stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            for it = subtrt'
                
%                 % input/output
%                 predictor = ([incentive,cost,instruction,effort_t,cost_t,trial_task,trial_hand,hand]);
%                 predictor = nanzscore(predictor(trt==treatmentList(it),:));
%                 indcat =  [1 2 3 5];
%                 predictor(:,indcat) = predictor(:,indcat).*(predictor(:,indcat)>0) ;
%                 y = [ effort(trt==treatmentList(it)) , rest(trt==treatmentList(it)) ];
%                 input = predictor';
%                 observation = y';
%                 
%                 
%                 % formula
%                 g_fname = @g_effortAccu;
% 
%                 % priors
%                 nphi = 10 ;
%                 param = struct;
%                 param.prior.mu =    ones(1,10);
%                 param.prior.sigma = ones(1,10);
%                 param.prior.sigma([1 8 9 10])=0; % fixed t0
% 
%                 param.type = repmat({'Phi'},1,nphi);
%                 param.labels = {'t0','tI','kE0','kED','kR0','kRI','kRD','kNt','kNh','kH'};
%                 param.transform.direct = [repmat({@identity},1,10)];
%                 inG.transform = param.transform.direct(ismember(param.type,'Phi'));
%                 opt.display=0;
%                 [priors]  = setParam(param,opt);
%                 
%                 
%                 % hyperpriors
%             
%                 % dimensions
%                 dim = struct('n',0,... % number of hidden states
%                         'n_theta',0 ,...   % number of evolution parameters
%                         'n_phi',nphi,...     % number of observation parameters
%                         'p',2,...          % output (data) dimension
%                         'n_t',size(observation,2));   % number of time samples or trials
%                 options.dim             = dim;
%     
%                 % options
%                 options.isYout  =  zeros(size(observation,1),size(observation,2));   % data exclusion
%                 observation(isnan(observation))=0;
%                 input(isnan(input))=0;
% 
%                 options.priors          = priors;   % include priors in options structure
%                 options.inG             = inG;      % input structure (grid)
%                 options.dim             = dim;
%                 options.DisplayWin      = 0;
%                 options.verbose         = 0;
%                 options.updateHP = 1 ;
%                 options.kernelSize = 1 ;
%                 options.extended=1;
%                 options.sources(1).out  = 1;    
%                 options.sources(1).type = 0;        % Normal continous data
%                 options.sources(2).out  = 2;    
%                 options.sources(2).type = 0;        % Normal continous data
% 
%                 % fit
%                 [posterior,inversion] = VBA_NLStateSpaceModel(observation,input,[],g_fname,dim,options);
%                 extract_model_output;
%                 disp(output.parameters);
%                 % extract
%                 Y(it,:,isub) = param.estimate;
                
                % select
                tab = result{isub}.gripAccu.inferential  ; 
                isess = find(tab.treatment==treatmentList(it));
                estimate = tab{:,4:13};
                % extract
                Y(it,:,isub) = estimate(isess,:);
                denominator = abs(Y(it,1,isub));
%                 Y(it,:,isub) = Y(it,:,isub)./denominator;
%                 Y(it,:,isub) = log(Y(it,:,isub));

                
            end
            
            % average normalization
%             Y(:,:,isub) = errorscore(Y(:,:,isub));

            denominator = abs(Y(2,:,isub));
            Y(:,:,isub) = Y(:,:,isub)./repmat(denominator,3,1);
%             denominator = abs(Y(2,1,isub));
%             Y(:,:,isub) = Y(:,:,isub)./denominator;
            Y2(:,:,isub) =  Y(:,:,isub)- repmat(Y(2,:,isub),3,1);
            

    end
    
    
    %% display glm
        fig = figure; set(fig,'Name','gripAccu_opioid');
            
            dsub=3;
%             y = nanmedian(Y,dsub);
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);
            
            hold on; clear h;
            for it=1:3
                x = [1:7];
%                 x = [1:nbin(2)];
                xx = x + (it-2)*0.2;
                [ h(it)  ] = barplot( xx ,y(it,x),z(it,x), col{it} );
% %                 [ h(it)  ] = barplot( x ,y(it),z(it), col{3} );
                h(it).BarWidth = 0.2;
                
%                 [ h ] = plotSpread(Y','distributionColors',col);
%                 m = findobj('-property','MarkerFaceColor');
%                 set(m,'MarkerSize',14);
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            effectName = {'t0','tI','kE0','kED','kR0','kRI','kRD','kNt','kNh','kH'};
%             effectName = {'intercept','incentive','effort_{t-1}','trial_{task}','trial_{hand}','hand','cost_{explicit}'};

            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = effectName(x);
            ylabel('parameters '); 
            title('model of effort allocation ');  
%             title('GLM of rest duration '); 
            ax.XLim = [0 x(end)+1];
            
    % format
    setFigProper('FontSize',20,'LineWidth',2);


%% learning
% valence effect

  nbin = [ 3 , 2 , nsub ];

    % variables
        Y = nan(nbin);
        Y2 = nan(nbin);


    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='learning') ;
            tab = tab(selection,:);
            
        % variables
%             correct = tab.isOptimalChoice;
            correct = tab.optimalChoice;
            money = tab.outcome*10;
            valence = tab.pairValence;

            trt = tab.treatment;
            trt = removecats(trt,'0');
            session = tab.sessionNumber;
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            trt = reordercats(trt,treatmentList(subtrt));
            ntrt = numel(subtrt);

            
            ysub = tools.tapply(correct,{trt,valence},@nanmean);   
%             ysub = tools.tapply(money,{trt,valence},@nansum);
            Y(subtrt,:,isub) =  ysub(:,[1 3]);
            Y2(subtrt,:,isub) =  Y(subtrt,:,isub)- repmat(Y(2,:,isub),ntrt,1);


    end
    
    % display
        fig = figure; set(fig,'Name','opioid_learning');
            
            dsub=3;
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);

            hold on;
            for it=1:3
                x = [0 5] + it;
                [ h(it)  ] = barplot( x ,y(it,:),z(it,:), col{it} );
%                 [ h(it)  ] = barplot( x ,y(it),z(it), col{3} );
                h(it).BarWidth = 0.2;
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            ax = gca; 
            ax.XTick = [2 7];
            ax.XTickLabel = {'loss','gain'};
            ylabel('correct choice (%) '); 
            ax.XLim = [0 x(end)+1];
            ax.YLim = [0.5 1];


        % format
        setFigProper('FontSize',20,'LineWidth',2);
        
%%
taskname = 'learning';
yname = 'trialDuration';

% ftransform = @ (x) normalize(x,'max');
ftransform = @identity;

fname = @nanmean;
% fname = @nanmax;

ylegend = 'force duration (% trial duration)';
stat = barcomp_motiscan_opioid(data,yname,ylegend,taskname,...
                               'fname',fname,'ftransform',ftransform,...
                               'sessionSplit',0); % without session interaction

%% learning
% learning effects

%   nbin = [ 3  , nsub ];
  nbin = [ 3 , 4 , nsub ];
%   nbin = [ 3 , 2 , nsub ];

    % variables
        Y = nan(nbin);
        Y2 = nan(nbin);


    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='learning') ;
            tab = tab(selection,:);
            stat = result{isub}.learning.inferential  ; 
            
        % variables
            correct = tab.isOptimalChoice;
            valence = tab.pairValence;
            
            BCA = stat.BCA;
%             ysub = [stat.repetition_previousLoss stat.repetition_previousNeutral_Loss ,...
%                     stat.repetition_previousNeutral_Gain , stat.repetition_previousGain ];
%             ysub = [ stat.dcorrect_Loss , stat.dcorrect_Gain ];
            ysub = [stat.correct_infirmed_Loss stat.correct_confirmed_Loss ,...
                    stat.correct_infirmed_Gain , stat.correct_confirmed_Gain ];
%             ysub = [ stat.beta_DV , stat.beta_DV ];
%             ysub = [ stat.beta_elapsedTime , stat.beta_elapsedTime ];
%             ysub = [ stat.repetition_Neutral , stat.repetition_Neutral ];
%             ysub = [ stat.correct_Loss , stat.correct_Gain ];

            trt = stat.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,treatmentList);


            
            session = tab.sessionNumber;
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            [~,indtrt] = sort(trt);
            ntrt = numel(subtrt);

%             ysub = tools.tapply(correct,{trt,valence},@nanmean);
%             Y(subtrt,:,isub) =  ysub(:,[1 3]);
%             Y2(subtrt,:,isub) =  Y(subtrt,:,isub)- repmat(Y(2,:,isub),ntrt,1);
            
            Y(subtrt,:,isub) =  ysub(indtrt,:);
            Y2(subtrt,:,isub) =  Y(indtrt,:,isub)- repmat(Y(indtrt(2),:,isub),ntrt,1);
            
            
            
            


    end
    
    % display
        fig = figure; set(fig,'Name','opioid_learning');
            
%             dsub=2;
            dsub=3;
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);

            hold on; clear h ; 
            for it=1:3
%                 x = 1;
                x = [1:nbin(2)];
                xx = x + 0.2*(it-2);
%                 [ h(it)  ] = barplot( xx ,y(it),z(it), col{it} );
                [ h(it)  ] = barplot( xx ,y(it,x),z(it,x), col{it} );
                h(it).BarWidth = 0.2;
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
%             varNames = {'BCA'};
%             varNames = {'repetition_previousLoss','repetition_previousNeutral_Loss',...
%                         'repetition_previousNeutral_Gain','repetition_previousGain'};
            varNames = {'correct infirmed Loss','correct confirmed Loss',...
                        'correct infirmed Gain','correct confirmed Gain'};
            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = varNames(x);
            ax.XTickLabelRotation = 45;
            ax.XLim = [0 x(end)+1];

        % format
        setFigProper('FontSize',20,'LineWidth',2);
        
%% learning effect

  nbin = [ 3 , 2 , 12, nsub ];

    % variables
        X = nan(nbin);
        Y = nan(nbin);
        Y2 = nan(nbin);

    for isub = 1:nsub
        try
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='learning') ;
            tab = tab(selection,:);
            
        % variables
%             correct = tab.isOptimalChoice;
            correct = tab.optimalChoice;
            prediction = tab.isOptimalChoice;
%             prediction = tab.predicted_isOptimalChoice;
            valence = tab.pairValence;
            block = tab.blockNumber + (tab.sessionNumber-1)*3;
            nt = tab.trialNumberByValence;
%             nt = tab.trialNumber;
            money = nan(size(nt));
            for i = unique(valence)'
               nt(valence==i) = quantileranks(nt(valence==i),nbin(3)) ;
               nt(nt==0) = NaN;
               for b = unique(block)'
                    money(valence==i & block==b) = cumsum(tab.outcome(valence==i & block==b)*10);
               end
            end

            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            session = tab.dayNumber;

            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            var = correct;
            xsub = tools.tapply(nt,{trt,valence,nt},@nanmean);
            ysub = tools.tapply(var,{trt,valence,nt},@nanmean);

            X(subtrt,:,:,isub) =  xsub(:,[1 3],:);
            Y(subtrt,:,:,isub) =  ysub(:,[1 3],:);

            ysub = tools.tapply(prediction,{trt,valence,nt},@nanmean);
            Y2(subtrt,:,:,isub) =  ysub(:,[1 3],:);

        end
    end
    
   %% display
        fig = figure; set(fig,'Name','gripAccu_learning2');
            
            dsub=4;
            x = nanmedian(X,dsub);
            y = nanmean(Y,dsub);
            y2 = nanmean(Y2,dsub);
            z = sem(Y,dsub);

            hold on;
            for it=1:3
                
                xx = permute(x(it,1,:),[3 1 2]);
                yy = permute(1-y(it,1,:),[3 1 2]);
                yy2 = permute(1-y2(it,1,:),[3 1 2]);
                zz = permute(z(it,1,:),[3 1 2]);
                [~,~,h] = errorscat(xx,yy,zz,col{it});
%                 h.LineStyle='none';
                h = plot(xx,yy2,'Color',col{it});
                h.LineStyle='--';
                
                xx = permute(x(it,2,:),[3 1 2]);
                yy = permute(y(it,2,:),[3 1 2]);
                yy2 = permute(y2(it,2,:),[3 1 2]);
                zz = permute(z(it,2,:),[3 1 2]);
                [~,~,h] = errorscat(xx,yy,zz,col{it});
%                 h.LineStyle='none';
                h = plot(xx,yy2,'Color',col{it});
            end
%             legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            ax = gca; 
            ax.XTick = xx;
            ax.XTickLabel = cellfun(@num2str,num2cell(2*xx),'UniformOutput',0);
            xlabel(' trial number '); 
            ylabel(' choice of best/worst option (gain/loss) '); 
            ax.XLim = [0 max(x(it,2,:))+1];
            ax.YLim = [0 1];
            
            
%% second-level parameters stats
  nbin = [ 3,4 , nsub ];

    % variables
        Y = nan(nbin);
        Y2 = nan(nbin);


    for isub = 1:nsub
        % select
            tab = result{isub}.learning.inferential  ; 
            
        % variables

            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            session = [1 2 3]';

            
        % stats
            [~,subtrt] = ismember(trt,treatmentList);
            ntrt = numel(subtrt);
            ysub = [tab.alpha tab.beta tab.kR tab.kP ];
%             ysub = errorscore([tab.alpha tab.beta tab.kR tab.kP ],2);
            Y(subtrt,:,isub) =  ysub;
            Y2(subtrt,:,isub) =  Y(subtrt,:,isub)- repmat(Y(2,:,isub),ntrt,1);
    end
    
    % display
        fig = figure; set(fig,'Name','learning_param');
            
            dsub=3;
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);

            hold on;
            for it=1:3
                x = [1 2 3 4];
                xx = x + (it-2)*0.2;
                [ h(it)  ] = barplot( xx ,y(it,:),z(it,:), col{it} );
                %                 [ h(it) ] = myboxplot( x , reshape(Y(it,:),nsub,1) , col{it} );
                h(it).BarWidth = 0.2;
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending$
            paramName = {'alpha','beta','kR','kP'};
            ax = gca; 
            ax.XTick = x;
            ax.XTickLabel = paramName(x);
            ylabel(' parametric effects (dev from the mean) '); 
            ax.XLim = [0 max(x)+1];
            
%% 3 ) Learning task
taskname = 'learning';

%% Data completion
for isub = 1:nsub % subject loop
        % select
        tab = data{isub}.battery.table; 
        selection =  (tab.task==taskname) ;
        tab = tab(selection,:);
        % variables
        accuracy = double(tab.isOptimalChoice == tab.optimalChoice);
        valence = tab.pairValence;
        money = tab.outcome*10;
        trt = tab.treatment;
        trt = removecats(trt,'0');
%         trt = reordercats(trt,treatmentList);
        % stats
        data{isub}.(taskname).table.accuracy = accuracy;
        [x,id,id2] = findgroups(valence(valence~=0),trt(valence~=0));
        y = accuracy(valence~=0);
        ycond = splitapply(@nanmean,y,x);
        result{isub}.(taskname).inferential.optimal_Gain = ycond(4:6);
        result{isub}.(taskname).inferential.optimal_Loss = ycond(1:3);
end

%% Parametric effect of treatments:

%   - predicted variables
% yname = {'correct_Gain','correct_Loss'};
yname = {'optimal_Gain','optimal_Loss'};
% yname = {'repetition_Neutral','repetition_previousLoss','repetition_previousNeutral_Loss','repetition_previousNeutral_Gain','repetition_previousGain'};
% yname = {'correct_infirmed_Loss','correct_confirmed_Loss','correct_infirmed_Gain','correct_confirmed_Gain'};
% yname = {'alpha','beta','kR','kP','BCA'};

%   - transform
fname = @nanmean;
% ftransform = @log;
ftransform = @identity;
% ftransform = @ (x) normalize(x,'variance');
%   - legend
ylegend = '';
% ylegend = 'parameters';

%   - order interaction
order = [];
%   - statistic:
%   'morphine_contrast','naloxone_contrast','opioid_contrast',...
%   'morphine_ttest','naloxone_ttest','opioid_ttest'
statistic = {'opioid_ttest'};
% statistic = {'naloxone_contrast','morphine_contrast'};
% statistic = {'naloxone_ttest','morphine_ttest'};
%   - plot: 'bar','dot','jitter','boxplot','violin','line'
plot = {'boxplot'};
sessionEffect=0;

[p,score,stat] = barcomp_stat_motiscan_opioid(result,yname,ylegend,taskname,...
                               'fname',fname,...
                               'ftransform',ftransform,...
                               'sessionSplit',sessionEffect,...
                               'order',order,...
                               'statisticalTest',statistic,...
                               'plotType',plot);
                           
%% Bivariate effect of treatment
%   - predicted variable
% yname = 'isOptimalChoice';
yname = 'accuracy';
% yname = 'RT';

%   - factor
xname = 'trialNumberByValence';
% xname = 'correctLogOddRatio';

interaction = 'pairValence';
%   - transform
ftransform = @identity;
%   - statistical function
fname = @nanmean;
% fname = @nanmax;
%   - session interaction
sessionEffect=0;
%   - legend
legend = {'trial (n)','correct choice(%)','valence'};
%   - order interaction
order = [];
% order = orderGrip;
%   - plot: 'dot','jitter','line'
% plot = {'dot','errorbar'};
% plot = {'jitter','errorbar'};
plot = {'errorbar'};

%   - interaction


figure;
bivarcomp_motiscan_opioid(data,xname,yname,legend,taskname,...
                               'fname',fname,...
                               'ftransform',ftransform,...
                               'sessionSplit',sessionEffect,...
                               'order',order,...
                               'statisticalTest',statistic,...
                               'plotType',plot,...
                               'interactionfactor',interaction);
            
            
            
%% volterra decomposition
% regressions
    tmax = 8; 
    valList = [-1 1];
    nbin = [ 3 , 2 , tmax , nsub ];

 % variables
    Y = nan(nbin);
    Y2 = nan(nbin);

    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='learning') ;
            tab = tab(selection,:);
            
        % variables
        
        for iv = [1 2]
            selection = (tab.pairValence==valList(iv));
            correct = tab.isOptimalChoice(selection);
            prediction = tab.predicted_isOptimalChoice(selection);
            valence = tab.pairValence(selection);
            outcome = tab.outcome(selection);
            nt = tab.trialNumberByValence(selection);
            
            trt = tab.treatment(selection);
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            session = tab.sessionNumber(selection);
            
            outcome_t = nan(numel(outcome),tmax);
            choice_t = nan(numel(outcome),tmax);            
            for t = 1:tmax
                outcome_t(t+1:end,t) = outcome(1:end-t);
                outcome_t(nt<=t,t) = NaN;
                choice_t(t+1:end,t) = (correct(1:end-t)*2)-1;
                choice_t(nt<=t,t) = NaN;
            end
            
            
            % first level stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            for it = subtrt'
                
%                 predictor = (outcome_t);
%                 predictor = (choice_t);
                predictor = (outcome_t.*choice_t);
                predictor = [(outcome_t.*choice_t) choice_t];
                y = correct;
                predictor = predictor(trt==treatmentList(it),:);  
                y = y(trt==treatmentList(it));
                
                formula = 'linear';
                stat = fitglm(predictor,y,formula,'Intercept',false);
                coef = stat.Coefficients;
                disp(coef);
               Y(it,iv,:,isub) = coef.Estimate(1:tmax);
               Y2(it,iv,:,isub) = coef.Estimate([1:tmax] + tmax );
            end
            
        end

    end
    
    
  % display
    fig = figure; set(fig,'Name','learning_volterra');

        dsub=4;
        var = Y2;
        y = nanmedian(var,dsub);
        z = sem(var,dsub);
        valence=2;
        y = reshape(y(:,valence,:),[3 tmax]);
        z = reshape(z(:,valence,:),[3 tmax]);

        hold on;
        for it=1:3
            x = [1:tmax];
            xx = x + (it-2)*0.2;
            [ h(it)  ] = barplot( xx ,y(it,:),z(it,:), col{it} );
            %                 [ h(it) ] = myboxplot( x , reshape(Y(it,:),nsub,1) , col{it} );
            h(it).BarWidth = 0.2;
        end
        legend([h(1) h(2) h(3)],treatmentList);

    % legending$
%         paramName = {'alpha','beta','kR','kP'};
        ax = gca; 
        ax.XTick = x;
%         ax.XTickLabel = paramName(x);
        xlabel('time lag : t_i'); 
        ylabel(' parametric effects of o(t-t_i)'); 
        title('volterra decomposition of choice(t)')
        ax.XLim = [0 max(x)+1];
            
%% rating task

  nbin = [ 3 , 6 , nsub ];

    % variables
        Y = nan(nbin);
        Y2 = nan(nbin);
        
        
    for isub = 1:nsub
        try
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='rating') ;
            tab = tab(selection,:);
            
        % variables
            rating = tab.rating./100;
            dim = tab.dimension;
            item = tab.itemNumber;
            type = tab.itemSubtype;
            
            trt = tab.treatment;
            trt = removecats(trt,'0');
            session = tab.sessionNumber;
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            trt = reordercats(trt,treatmentList(subtrt));
            ntrt = numel(subtrt);
            
%             ysub = tools.tapply(rating,{trt,dim},@nanmean);
            ysub = tools.tapply(rating,{trt,dim,type},@nanmean);

            ysub = tools.tapply(rating,{trt,dim},@nanvar);


            Y(subtrt,:,isub) =  ysub;
            Y2(subtrt,:,isub) =  Y(subtrt,:,isub)- repmat(Y(2,:,isub),ntrt,1);
        end
    end
    
    % display
        fig = figure; set(fig,'Name','opioid_rating');
            
            dsub=3;
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);

            for idim=1:6
                subplot(2,3,idim);
                hold on; 
                h = gobjects(1,3);
                for it=1:3
                    x = it;
                    [ h(it)  ] = barplot( x ,y(it,idim,1),z(it,idim,1), col{it} );
%                     [ h(it) ] = myboxplot( x , reshape(Y(it,idim,:),nsub,1) , col{it} );
                    h(it).BarWidth = 0.5;
                end
                
                legend([h(1) h(2) h(3)],treatmentList);

            % legending
                ax = gca; 
                ax.XTick = [];
                ylabel(' average ratings  '); 
%                 ylabel(' ratings variance  '); 
                ax.XLim = [0 max(x)+1];
%                 ax.YLim = [0 1];
                ax.YLim = [0 1];

                title(char(dimensionList(idim + (idim>=2) )))
            end
            
            
        % format
        setFigProper('FontSize',20,'LineWidth',2);
            
%% weight task

%   nbin = [ 3 , 3 , nsub ];
%   nbin = [ 3 , 7 , nsub ];
  nbin = [ 3 , 14 , nsub ];

    % variables
        Y = nan(nbin);

    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='rating') ;
            tab = tab(selection,:);
            
        % variables
            rating = tab.rating./100;
            dim = tab.dimension;
            item = tab.itemNumber;
                        
            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            session = tab.sessionNumber;
        
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            ntrt = numel(subtrt);
            ysub = tools.tapply(rating,{trt,dim,item},@nanmean);
            gainRating = reshape(ysub(:,4,1:7),ntrt,7);
            lossRating = - reshape(ysub(:,5,1:7),ntrt,7);
            ysub = [fliplr(lossRating) , gainRating];
            Y(subtrt,:,isub) =  ysub;            
%             Y(subtrt,:,isub) =  errorscore(ysub);

    end
    
    % display
        fig = figure; set(fig,'Name','opioid_rating_money');
            
            dsub=3;
            levels =  [-20 -10 -5 -1 -0.50 -0.20 -0.01 0.01 0.20 0.50 1 5  10 20];
            y = nanmean(Y,dsub);
            z = sem(Y,dsub);

            for idim=1:3
                hold on; 
                h = gobjects(1,3);
                for it=1:3
                    x = [1:nbin(2)];
                    xx = (x)-7.5;
                    [~,~,h(it)] = errorscat(xx ,y(it,:,:),z(it,:,:), col{it});
                end
                legend([h(1) h(2) h(3)],treatmentList);

            % legending
                ax = gca; 
                ax.XScale='linear';
                ax.XTick = xx ;
                ax.XTickLabel = cellfun(@num2str,num2cell(levels),'UniformOutput',0) ;

                xlabel(' incentive (€) '); 
                ylabel(' rating (%) '); 
%                 ax.XLim = [0 max(x)+1];
%                 ax.YLim = [0 1];
            end
            
    % format
    setFigProper('FontSize',20,'LineWidth',2);

            
%% weight task

  nbin = [ 3 , 3 , nsub ];

    % variables
        Y = nan(nbin);

    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='weight') ;
            tab = tab(selection,:);
            
        % variables
            accept = tab.isGoChoice;
            dim = tab.dimension;
            
            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            session = tab.sessionNumber;

            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            ysub = tools.tapply(accept,{trt,dim},@nanmean);
            Y(subtrt,:,isub) =  ysub;            

    end
    
    % display
        fig = figure; set(fig,'Name','opioid_weight');
            
            dsub=3;
            y = nanmedian(Y,dsub);
            z = sem(Y,dsub);

            for idim=1:3
                subplot(1,3,idim);
                hold on; 
                h = gobjects(1,3);
                for it=1:3
                    x = it;
                    [ h(it)  ] = barplot( x ,y(it,idim,1),z(it,idim,1), col{it} );
                     h(it).BarWidth = 0.5;
%                     [ h(it) ] = myboxplot( x , reshape(Y(it,idim,:),nsub,1) , col{it} );
                end
                legend([h(1) h(2) h(3)],treatmentList);

            % legending
                ax = gca; 
                ax.XTick = [];
                ylabel(' acceptance rate (%) '); 
                ax.XLim = [0 max(x)+1];
                ax.YLim = [0 1];
                title(char(dimensionList(idim + 7 )))
            end

            
  % second level stats
    accept = reshape(Y,nsub*3*3,1);
    trt = repmat([1 2 3]',nsub*3,1);
    dim = repmat([1 1 1 2 2 2 3 3 3]',nsub,1);
    trt = nominal(treatmentList(trt))';
    trt = addlevels(trt,{'0'}); trt = reorderlevels(trt,{'0',treatmentList{:}});
    dim = nominal(dimensionList(dim+7))';
    dim = addlevels(dim,{'0'}); dim = reorderlevels(dim,['0',cellstr(dimensionList(8:10))]);

    % glm
    formula = 'accept ~ -1 + dim + dim:trt';
    stat = fitglm(table(dim,trt,accept),...
                    formula,...
                   'CategoricalVars',[1 2]);
    coef = stat.Coefficients;
    disp(coef);
    
%     design = [ double(dim),double(trt) ];
%     design(isnan(design))=0;
%     accept(isnan(accept))=0;
%     [p,stat] = GLM_contrast([design ,design(1).*design(2)],accept,[0 0 1]','F',0);
    
    % anovan
%     [p,~,stat] = anovan(accept,{dim trt},...
%                     'model','interaction');
               
    [p,F,d] = coefTest(stat,[0 0 0  0 0 0  0 0 0  1 0 1 ] );
                         
%% acceptance function
nbin = [ 3 , 3 , 6 , nsub ];

    % variables
        X = nan(nbin);
        Y = nan(nbin);

    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='weight') ;
            tab = tab(selection,:);
            
        % variables
            accept = tab.isGoChoice;
            dim = tab.dimension;
%             value = tab.benefitItemRatingValue./100;
            value = tab.costItemRatingValue./100;
            qvalue = value;
            for i = unique(dim)'
               qvalue(dim==i) = quantileranks(value(dim==i),nbin(3)) ;
               qvalue(qvalue==0) = NaN;
            end


            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            session = tab.sessionNumber;

            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            xsub = tools.tapply(value,{trt,dim,qvalue},@nanmean);
            ysub = tools.tapply(accept,{trt,dim,qvalue},@nanmean);
            X(subtrt,:,:,isub) =  xsub;            
            Y(subtrt,:,:,isub) =  ysub;            

    end
    
    % display
        fig = figure; set(fig,'Name','opioid_weight');
            
            dsub=4;
            x = nanmean(X,dsub);
            y = nanmean(Y,dsub);
            z = sem(Y,dsub);

            for idim=1:3
                subplot(1,3,idim);
                hold on; 
                h = gobjects(1,3);
                for it=1:3
                    xx = reshape(x(it,idim,:,1) , nbin(3) , 1);
                    yy = reshape(y(it,idim,:,1) , nbin(3) , 1);
                    zz = reshape(z(it,idim,:,1) , nbin(3) , 1);
                    [~,~,h(it)] = errorscat(xx,yy,zz,col{it});
                end
                legend([h(1) h(2) h(3)],treatmentList);

            % legending
                ax = gca; 
%                 xlabel(' benefit value '); 
                xlabel(' cost value '); 
                ylabel(' acceptance rate (%) '); 
                ax.XLim = [0 1];
                ax.YLim = [0 1];
                title(char(dimensionList(idim + 7 )))
            end
            
%% second-level parameters stats
  nbin = [ 3 , 3 , nsub ];

    % variables
        Y = nan(nbin);

    for isub = 1:nsub
        try
        % select
            tab = result{isub}.battery.inferential  ; 
            
        % variables

            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            session = [1 2 3]';

            
        % stats
            [~,subtrt] = ismember(trt,treatmentList);
%             ysub = [tab.muR  tab.muE  tab.muE];
%             ysub = [tab.stdR  tab.stdP  tab.stdE];
            ysub = [tab.muR./tab.stdR  tab.muE./tab.stdP  tab.muE./tab.stdE];
%             ysub = [tab.alpha  tab.bR  tab.theta];

            Y(subtrt,:,isub) =  ysub;
        end

    end
    
    % display
        fig = figure; set(fig,'Name','pref_param');
            
            dsub=3;
            y = nanmedian(Y,dsub);
            z = sem(Y,dsub);
            ylist = {'muR','muP','muE'};
            
            for idim=1:3
                subplot(1,3,idim);
                hold on;
                for it=1:3
                    x = it;
                    [ h(it)  ] = barplot( x ,y(it,idim),z(it,idim), col{it} );
    %                 [ h(it) ] = myboxplot( x , reshape(Y(it,:),nsub,1) , col{it} );
                    h(it).BarWidth = 0.5;
                end
                
                legend([h(1) h(2) h(3)],treatmentList);            
                % legending
                ax = gca; 
                ax.XTick = [];
                title(ylist{idim}); 
            end
            
%% 4 ) Preference tasks
taskname = 'battery';


%% Bivariate effect of treatment
taskname = 'rating';
% taskname = 'weight';


%   - predicted variable
yname = 'rating';
% yname = 'isGoChoice';

%   - factor
% xname = 'dimension';
% xname = 'benefitItemRatingValue';
xname = 'itemSubtype';

interaction = [];
% interaction = 'dimension';
% interaction = 'pairValence';
%   - transform
ftransform = @identity;
%   - statistical function
fname = @nanmean;
% fname = @nanmax;
%   - session interaction
sessionEffect=0;
%   - legend
legend = {'dimension','rating(%)'};
% legend = {'benefit rating','acceptance(%)','dimensions'};
statistic = {'opioid_ttest'};
%   - order interaction
order = [];
% order = orderGrip;
%   - plot: 'dot','jitter','line'
% plot = {'dot','errorbar'};
% plot = {'jitter','errorbar','boxplot'};
plot = {'boxplot'};
% plot = {'errorbar'};

%   - interaction


figure;
bivarcomp_motiscan_opioid(data,xname,yname,legend,taskname,...
                               'fname',fname,...
                               'ftransform',ftransform,...
                               'sessionSplit',sessionEffect,...
                               'order',order,...
                               'statisticalTest',statistic,...
                               'plotType',plot,...
                               'interactionfactor',interaction);
                           
%% Effect on item-specific ratings


% --- 1. Item-specific additive treatment effect
% variables
ndim=3;
ntrt=3;
nitem=24;
DIM = [];
T = [];
S = [];
Y = nan(ndim,ntrt,nsub);
Y2 = nan(ndim,ntrt,nsub);

for isub = 1:nsub
    % select
        tab = data{isub}.rating.table  ; 
        subset = (tab.dimension=='writtenReward' |tab.dimension=='writtenPunishment' |tab.dimension=='writtenEffort' );
        tab = tab( subset,:);

    % variables
        trt = tab.treatment;
        trt = removecats(trt,'0');
        trt = reordercats(trt,treatmentList);
        rating = tab.rating./100;
        dim = tab.dimension;
        dim = removecats(dim,dimensionList([2 5:14]));
        dim = reordercats(dim,dimensionList([1 3 4]));
        item = tab.itemNumber;
        sess = tab.sessionNumber;

    % first level stats
        [g,ind_dim,ind_item,ind_trt] = findgroups(dim,item,trt);
        ysub = splitapply(@unique,rating,g);
        
        ysub = permute(reshape(ysub,[3 24 3]),[3 2 1]);
        ysub2 = ysub - mean(ysub,3);
        ysub = nanmean(ysub,2);
        ysub2 = nanmean(ysub2,2);
        Y(:,:,isub) =   reshape(ysub,[ndim,ntrt]);
        Y2(:,:,isub) =   reshape(ysub2,[ndim,ntrt]);
        
        [g,ind_trt,ind_dim] = findgroups(trt,dim);
        subsess = splitapply(@nanmean,sess,g);
        S = [S ; subsess];
        DIM = [DIM ; ind_dim];
        T = [T ; ind_trt];

end
    
% second level stats
% reshape
y = reshape(Y2,ndim*ntrt*nsub,1);
dim = reshape(DIM,ndim*ntrt*nsub,1);
trt = reshape(T,ndim*ntrt*nsub,1);
sess = reshape(S,ndim*ntrt*nsub,1);

% glm
% predictor = [ (trt=='naloxone') , (trt=='placebo'), (trt=='morphine') , dim, sess ];
% formula = ['y ~ -1 + placebo + naloxone + morphine + dim:(placebo + naloxone + morphine) + session '];
predictor =  table((T=='naloxone')*(-1) + (T=='placebo')*0 + (T=='morphine')*1 , dim, sess );
predictor.Properties.VariableNames = {'opioid','dim','session'};
formula = ['y ~ -1 + opioid + dim:opioid + session'];
stat = fitglm([predictor,table(y)],formula);
coef = stat.Coefficients;
disp(coef);


% display
clear g;
alpha=0.8;
g = gramm('x',dim,'y',y,'color',trt);
g.set_color_options('map',vertcat(col{:}),'lightness',100);
g.set_order_options('x',dimensionList([1 3 4]),'color',treatmentList);
% g.geom_jitter('height',0.01);
g.stat_boxplot('width',0.9);
g.set_names('x','','y','rating deviation from mean (%)','color','treatment');
g.axe_property('YLim',[mean(min(y)) mean(max(y))] + [-0.2 0.2]*mean(mean(y)));
% g.axe_property('XTickLabelRotation',45);
g.draw;
axes(g.facet_axes_handles);
set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                    'LineWidth',1.5,'FaceAlpha',alpha);

%% --- 2. Item-specific multiplicative treatment effect
% variables
ndim=3;
ntrt=3;
nitem=24;
DIM = [];
T = [];
S = [];
ITEM = [];
Y = nan(ndim,nitem,ntrt,nsub);
Y2 = nan(ndim,nitem,ntrt,nsub);

for isub = 1:nsub
    % select
        tab = data{isub}.rating.table  ; 
        subset = (tab.dimension=='writtenReward' |tab.dimension=='writtenPunishment' |tab.dimension=='writtenEffort' );
        tab = tab( subset,:);

    % variables
        trt = tab.treatment;
        trt = removecats(trt,'0');
        trt = reordercats(trt,treatmentList);
        rating = tab.rating./100;
        dim = tab.dimension;
        dim = removecats(dim,dimensionList([2 5:14]));
        dim = reordercats(dim,dimensionList([1 3 4]));
        item = tab.itemNumber;
        sess = tab.sessionNumber;

    % first level stats
        [g,g_trt,g_item,g_dim] = findgroups(trt,item,dim);
        ysub = splitapply(@unique,rating,g);        
        ysub = reshape(ysub,[ndim,nitem,ntrt]);
        ysub2 = ysub - mean(ysub,3);
        Y(:,:,:,isub) =   repmat(mean(ysub,3),1,1,3);
        Y2(:,:,:,isub) =   ysub2;
        
        subsess = splitapply(@nanmean,sess,g);
        S = [S ; subsess];
        DIM = [DIM ; g_dim];
        T = [T ; g_trt];
        ITEM = [ITEM ; g_item];

end
    
% second level stats
YGROUP = mean(mean(Y,4),3);
Y3 = Y2 + YGROUP;

% reshape
x = reshape(Y,ndim*ntrt*nitem*nsub,1);
y = reshape(Y3,ndim*ntrt*nitem*nsub,1);
dim = reshape(DIM,ndim*ntrt*nitem*nsub,1);
item = reshape(ITEM,ndim*ntrt*nitem*nsub,1);
trt = reshape(T,ndim*ntrt*nitem*nsub,1);
sess = reshape(S,ndim*ntrt*nitem*nsub,1);
itemRatingsByTrt = Y;
itemRatings = squeeze(mean(Y,3));

% display
dimList = dimensionList([1 3 4]);
titleList = {'reward','punishment','effort'};
clear g;
alpha=0.8;
for idim=1:3
    [~,rank] = sort(YGROUP(idim,:));
    g(1,idim) = gramm('x',nominal(item),'y',y,'color',trt,'subset',(dim==dimList{idim}));
    g(1,idim).set_color_options('map',vertcat(col{:}),'lightness',100);
    g(1,idim).set_order_options('x',nominal(rank),'color',treatmentList);
    g(1,idim).stat_summary('type','sem','geom',{'point','errorbar','line'});
%     g(1,idim).set_names('x','item number','y','rating average | treatment, item (%)','color','treatment');
    g(1,idim).set_names('x','item rank (sorted by average rating)','y','rating average | treatment, item (%)','color','treatment');
    g(1,idim).set_title(titleList{idim});
    g(1,idim).axe_property('XLim',[0 25],'YLim',[0 1],'XTick',[]);
    
end
g.draw;
for idim=1:3
axes(g(1,idim).facet_axes_handles);
set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                    'LineWidth',1.5,'FaceAlpha',alpha);
end

                           
%%

% taskname = 'battery';
taskname = 'weight';

% % model-based
% yname = {'muR','muP','muE','alpha','bR','t0','theta'};
% % ftransform = @log;
% ftransform = @identity;
% % ftransform = @ (x) normalize(x,'variance');
% 
% [p,F,stat] = barcomp_stat_motiscan_opioid(result,yname,'parameters',taskname,...
%                                'fname',@nanmean,'ftransform',ftransform,...
%                                'sessionSplit',0); % without session interaction

% model-based
yname = {'acceptanceRate_RE','acceptanceRate_PE','acceptanceRate_RP'};
%   - transform
fname = @nanmean;
% ftransform = @log;
ftransform = @identity;
% ftransform = @ (x) normalize(x,'variance');
%   - legend
ylegend = '';
% ylegend = 'parameters';

%   - order interaction
order = [];
%   - statistic:
%   'morphine_contrast','naloxone_contrast','opioid_contrast',...
%   'morphine_ttest','naloxone_ttest','opioid_ttest'
statistic = {'opioid_ttest'};
% statistic = {'naloxone_contrast','morphine_contrast'};
% statistic = {'naloxone_ttest','morphine_ttest'};
%   - plot: 'bar','dot','jitter','boxplot','violin','line'
plot = {'boxplot'};
sessionEffect=0;

[p,score,stat] = barcomp_stat_motiscan_opioid(result,yname,ylegend,taskname,...
                               'fname',fname,...
                               'ftransform',ftransform,...
                               'sessionSplit',sessionEffect,...
                               'order',order,...
                               'statisticalTest',statistic,...
                               'plotType',plot);
                           
%% logistic regression effect of opioids

% variables
nparam=6;
ntrt=3;
T = [];
S = [];
Y = nan(nparam,ntrt,nsub);
Y2 = nan(nparam,ntrt,nsub);

for isub = 1:nsub
    % select
    tab = data{isub}.weight.table  ; 

    % variables
    trt = tab.treatment;
    trt = removecats(trt,'0');
    trt = reordercats(trt,treatmentList);
    dim = tab.dimension;
    dim = removecats(dim,dimensionList([1:7,11:14]));
    dim = reordercats(dim,dimensionList([8 9 10]));
    % rating values|trt
%         reward = tab.rewardRatingValue./100;
%         punishment = tab.punishmentRatingValue./100;
%         punishment(dim=='PunishmentEffort') = - punishment(dim=='PunishmentEffort');
%         effort = tab.effortRatingValue./100;
    % rating values (average)
        benefit_item = tab.benefitItemNumber;
        cost_item = tab.costItemNumber;
        reward = zeros(size(trt));
        reward(dim=='RewardEffort') = itemRatings(1,benefit_item(dim=='RewardEffort'),isub);
        reward(dim=='RewardPunishment') = itemRatings(1,benefit_item(dim=='RewardPunishment'),isub);
        punishment = zeros(size(trt));
        punishment(dim=='PunishmentEffort') = - itemRatings(2,benefit_item(dim=='PunishmentEffort'),isub);
        punishment(dim=='RewardPunishment') = itemRatings(2,cost_item(dim=='RewardPunishment'),isub);
        effort = zeros(size(trt));
        effort(dim=='RewardEffort') = itemRatings(3,cost_item(dim=='RewardEffort'),isub);
        effort(dim=='PunishmentEffort') = itemRatings(3,cost_item(dim=='PunishmentEffort'),isub);
    sess = tab.sessionNumber;
    accept = tab.isGoChoice;

    % first level stats
    predictor = table(trt,dim,reward,punishment,effort,sess,accept);
    formula = ['accept ~ dim + reward + punishment + effort'];
    for i = 1:numel(treatmentList)
        subset = (trt==treatmentList{i});
        subpredictor = predictor(subset,:);
        stat = fitglm(subpredictor,formula);
        coef = stat.Coefficients;
        coef.Properties.RowNames{1} = 'dim_RewardEffort';
        Y(:,i,isub) = coef.Estimate;
    end
    Y2(:,:,isub) = Y(:,:,isub) - mean(Y(:,:,isub),2);
    g = findgroups(trt);
    subsess = splitapply(@unique,sess,g);
    S = [S,subsess];    

end
    
% second level stats
% reformat
Y2 = Y2 + nanmean(nanmean(Y,3),2);
paramNames = coef.Properties.RowNames;
y = reshape(Y2,nparam,ntrt*nsub,1)';
trt = repmat(nominal(treatmentList)',nsub,1);
sess = reshape(S,ntrt*nsub,1);

% glm
predictor = table((trt=='naloxone'),(trt=='placebo'),(trt=='morphine'),sess);
predictor.Properties.VariableNames = {'naloxone','placebo','morphine','session'};
formula = ['y ~ -1 + placebo + naloxone + morphine + session '];
contrast = [-1 0 1 0];
% predictor =  table((trt=='naloxone')*(-1) + (trt=='placebo')*0 + (trt=='morphine')*1, sess );
% predictor.Properties.VariableNames = {'opioid','session'};
% formula = ['y ~ 1 + opioid + session'];
for ip = 1:nparam
    stat = fitglm([predictor,table(y(:,ip),'VariableNames',{'y'})],formula);
    coef = stat.Coefficients;
    disp(paramNames{ip});
    disp(coef);
    [p,score,df] = coefTest(stat,contrast);
    fprintf('contrast: morphine > naloxone , F = %d, p = %d \n',score,p);
    [p,t] = linhyptest(stat.Coefficients.Estimate,stat.CoefficientCovariance,0,[-1 0 1 0],stat.DFE);
    fprintf('contrast: morphine > naloxone , F = %d, p = %d \n',t,p);
end

% reshape
y = reshape(y,nparam*ntrt*nsub,1);
param = reshape(repmat(nominal(paramNames'),ntrt*nsub,1),nparam*ntrt*nsub,1);
trt = repmat(trt,nparam,1);
sess = repmat(sess,nparam,1);
sub = reshape(repmat([1:nsub],nparam*ntrt,1),nparam*ntrt*nsub,1);

% display
clear g;
alpha=0.8;
g = gramm('x',param,'y',y,'color',trt);
g.set_color_options('map',vertcat(col{:}),'lightness',100);
g.set_order_options('color',treatmentList);
% g.stat_summary('type','sem','geom',{'bar','black_errorbar'});
g.stat_boxplot('width',0.9);
g.set_names('x','','y','logistic regression coefficients (au.)','color','treatment');
g.axe_property('XTickLabelRotation',45);
g.draw;
axes(g.facet_axes_handles);
set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                    'LineWidth',1.5,'FaceAlpha',alpha);
