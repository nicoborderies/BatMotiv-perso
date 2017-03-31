%% analyze_batmotiv_opioid_2level

% reset
clc;
clear all;
close all;

%% Prepare data analysis
%-------------------------------

% load dataset
datadir = 'B:\nicolas.borderies\projets\batmotiv\resultats\OPIOID';
codedir = 'B:\nicolas.borderies\projets\batmotiv\code.perso';
cd(datadir);
analysisname = ['motiscan_opioid_dataset.mat'];
load(analysisname);
cd(codedir);


        
%% Data analysis
%-------------------------------

% define analysis parameters
data = [groupData.OPIOID.subject ];
result = [groupResult.OPIOID.subject ];
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

%% 1) grip task
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
            xlabel('incentive(�)'); 
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

% parameters
nbin = [ ntrt , nsub ];


% Effect of Treatments:
%   - effort duration

    % prepare variables
        Y = nan(nbin);
        Y2 = nan(nbin);
    for isub = 1:nsub % subject loop
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='gripAccu') ;
            tab = tab(selection,:);
            
        % variables
            effort = tab.effortDuration;
%             effort = effort./max(effort(tab.treatment=='placebo'));
%             effort = normalize(effort,'zscore');
%             effort = tab.predicted_effortDuration;
            rest = tab.restDuration;
            totaltime = tab.effortDuration + tab.restDuration ;
            gain = tab.gain;
            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,treatmentList);
            session = tab.sessionNumber;
            
        % stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            trt = reordercats(trt,treatmentList(subtrt));

            
%             ysub = tools.tapply(effort,{trt},@nanmean);
%             ysub = tools.tapply(rest,{trt},@nanmean);
            ysub =  tools.tapply(effort,{session},@nanmean);
%             ysub =  tools.tapply(gain,{trt},@nanmax);

            Y(subtrt,isub) =  ysub;
            Y2(subtrt,isub) =  Y(subtrt,isub)- Y(2,isub);
     end
    
    
    % display
        fig = figure; set(fig,'Name','gripAccu_opioid');
            
            dsub=2;
            y = nanmean(Y,dsub);
            z = sem(Y2,dsub);
            
            hold on; clear h;
            for it=1:3
                x = it;
                [ h(it)  ] = barplot( x ,y(it),z(it), col{it} );
% %                 [ h(it)  ] = barplot( x ,y(it),z(it), col{3} );
                h(it).BarWidth = 0.5;
                
%                 [ h ] = plotSpread(Y','distributionColors',col);
%                 m = findobj('-property','MarkerFaceColor');
%                 set(m,'MarkerSize',14);
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
        % legending
            ax = gca; 
            ax.XTick = [];
            ylabel('force duration (sec) '); 
%             ylabel('rest duration (sec) '); 
%             ylabel('total gain (�)');
            ax.XLim = [0 x(end)+1];
            ax.YLim = [min(y) max(y)] + [-0.2 0.2]*mean(y);

        % format
        setFigProper('FontSize',20,'LineWidth',2);


%% cost effect
    nbin = [ 3 , 2 , nsub ];

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
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='gripAccu') ;
            tab = tab(selection,:);
            
        % variables

%             effort = tab.effortDuration;
            effort = tab.predicted_effortDuration;
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
            xlabel('incentive level (�)  '); 
            ylabel('force duration (sec) '); 
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
            xlabel('effort cummulated_{t-1}  (%fmax.sec) '); 
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
% raw effect
    nbin = [ 3 , 2 , nsub ];

 % variables
    Y = nan(nbin);

    for isub = 1:nsub
                
        % select
            tab = data{isub}.battery.table; 
            selection =  (tab.task=='gripAccu') ;
            tab = tab(selection,:);
            
        % variables
            effort = (tab.effortDuration);
%             effort = errorscore(tab.effortDuration);
%             effort = errorscore(tab.restDuration);

            cost = tab.costValue;
            incentive = tab.incentiveValue;

            trt = tab.treatment;
            trt = removecats(trt,'0');
            trt = reordercats(trt,sort(categories(trt)));
            
            session = tab.sessionNumber;
            
        % first level stats
            [~,subtrt] = ismember(unique(trt),treatmentList);
            Y(subtrt,1,isub) =  tools.tapply(effort,{trt},@nanmean);
            Y(subtrt,2,isub) =  tools.tapply(session,{trt},@unique);


    end
    
    % second level stats
%     trt = repmat([1 2 3]',nsub,1);
    trt = repmat([1 -1 0]',nsub,1);
    session = reshape(Y(:,2,:),nsub*3,1);
    effort = reshape(Y(:,1,:),nsub*3,1);
    
    formula = 'effort ~ 1 + trt + session';
    stat = fitglm([trt,session],effort,...
                   formula,...
                   'VarNames',{'trt','session','effort'},...
                   'CategoricalVars',1);
    coef = stat.Coefficients;
    coef.Properties.RowNames = {'morphine','naloxone','placebo','session'};
               
%     formula = 'effort ~ -1 + trt + session';
%     stat = fitglm([trt,session],effort,...
%                    formula,...
%                    'VarNames',{'trt','session','effort'});
%     coef = stat.Coefficients;
%     coef.Properties.RowNames = {'treatment','session'};

    disp(coef);

%     writetable(coef,'stat_gripAccu.xlsx','WriteRowNames',1);

    [p,F,d] = coefTest(stat,[1 1 1 0])
%     [p,F,d] = coefTest(stat,[0 0 0 1])
    [p,F,d] = coefTest(stat,[1 0 -1 0])


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
            correct = tab.isOptimalChoice;
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
            correct = tab.isOptimalChoice;
            prediction = tab.predicted_isOptimalChoice;
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
                [~,~,h(it)] = errorscat(xx,yy,zz,col{it});
                h(it).LineStyle='none';
                h(it) = plot(xx,yy2,'Color',col{it});
                h(it).LineStyle='--';
                
                xx = permute(x(it,2,:),[3 1 2]);
                yy = permute(y(it,2,:),[3 1 2]);
                yy2 = permute(y2(it,2,:),[3 1 2]);
                zz = permute(z(it,2,:),[3 1 2]);
                [~,~,h(it)] = errorscat(xx,yy,zz,col{it});
                h(it).LineStyle='none';
                h(it) = plot(xx,yy2,'Color',col{it});
            end
            legend([h(1) h(2) h(3)],treatmentList);
            
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

                xlabel(' incentive (�) '); 
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
