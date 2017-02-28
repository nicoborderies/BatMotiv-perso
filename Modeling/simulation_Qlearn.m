%% Simulation: Q-Value Reinforcement Learning Model

clear all;
clc;
close all;
 
 
 


% --- VBA model simulation --- %

% 1/ simulation environment
    % empirical parameters
    
    blockNumber = 2;
    trialNumber = 48;
    n_t = blockNumber*trialNumber;  % number of observations  

%     % experimental factor
%     block = NaN(1,trialNumber*blockNumber);
%     for iBlock = 1:blockNumber
%         block(1 + trialNumber*(iBlock-1):trialNumber*(iBlock))= iBlock*ones(1,trialNumber);
%     end
%     trial = repmat([1:trialNumber],1,blockNumber);
%     pairValence  = repmat([randperm(trialNumber/2) randperm(trialNumber/2)],1,blockNumber); % valence state    
%     lottery=ones(1,n_t);
%     for i=[1:3,(trialNumber/4+1:trialNumber/4+3)]
%         lottery(pairValence==i)=-1;                                             % -1 = unlikely outcome (25%), % 1 = likely outcome (75%)
%     end
%     pairValence = -2*mod(pairValence-1,2)+1;
%             
%     % initial conditions
%     initialTrialLoss =[];initialTrialGain=[];finalTrialLoss=[];finalTrialGain=[];
%     for iBlock = 1:blockNumber
%         initialTrialLoss(iBlock) = find(pairValence(block==iBlock) == -1,1,'first') + (iBlock-1)*trialNumber;
%         initialTrialGain(iBlock) = find(pairValence(block==iBlock) == 1,1,'first') + (iBlock-1)*trialNumber;
%         finalTrialLoss(iBlock) = find(pairValence(block==iBlock) == -1,1,'last') + (iBlock-1)*trialNumber;
%         finalTrialGain(iBlock) = find(pairValence(block==iBlock) == 1,1,'last') + (iBlock-1)*trialNumber;
%     end
% 
%     previousValence = zeros(1,trialNumber*blockNumber);
%     previousValence(setdiff(find(trial),min(initialTrialLoss,initialTrialGain)))  = pairValence( ~ismember(find(trial),max(finalTrialLoss,finalTrialGain)));
%     previousChoice = zeros(1,trialNumber*blockNumber);
%     previousOutcome = zeros(1,trialNumber*blockNumber);
% 
%     bin=zeros(size(trial));
%     nTrialByValence=numel(trial)/(numel(unique(block))*2);
%     nBin=6;
%     binSize = round(nTrialByValence/nBin);
%     nBin=nTrialByValence/binSize;
%     bin(pairValence==-1)=repmat(sort(repmat(1:nBin, 1, binSize)),1, numel(unique(block)));
%     bin(pairValence==1)=repmat(sort(repmat(1:nBin, 1, binSize)),1, numel(unique(block)));
% 
%     u = [  pairValence ; previousValence ; previousChoice ; previousOutcome ];                    
%     designList = {'','',''};
%     designMatrix{iD} = [x1 ; x2 ; x3];

    % model space 
%     modelList = {'','',''};
    f_fname = @f_Qlearn; % evolution function (Q-value learning)
    g_fname = @g_softmax; % observation function (softmax selection)
    h_fname = @h_stochasticOutcome; % feedback function (stochastic 2x2 state-action gain/loss space)

    % specifications    
%     specificationList = {'','',''};
%     fb.inH.lottery = lottery;
%     fb.inH.pairValence = pairValence;
%     fb.h_fname = h_fname;        
%     fb.indy = 3;
%     fb.indfb = 4;

    % model options
    inF.counterfactual = 0;
    inF.stateDependantLR = 0;
    inF.outcomeWeight = 2;
    inG.stateDependantT = 0;

    % model dimensions
        dim = struct('n',4,... % number of hidden states
            'n_theta',1 + inF.stateDependantLR + inF.outcomeWeight ,...    % number of evolution parameters
            'n_phi',1+inG.stateDependantT,...                             % number of observation parameters
            'p',1,...          % output (data) dimension
            'n_t',n_t);   % number of time samples or trials
    
    % parameters
        % model parameters
        paramList = {'alpha','',''};
            % Observation parameters 
            inverse_temperature = [3 ; 3];
%             inverse_temperature = [2 ; 4];

            inverse_temperature = log(inverse_temperature);
            for iP=1:2
                phi{iP} = inverse_temperature(iP);
            end
            priors.muPhi = phi{iP};
            priors.SigmaPhi = 0*eye(dim.n_phi);

            % Evolution parameters 
            alpha = [0.10 ; 0.50 ] ;
%             alpha = [0.35 ; 0.35 ] ;

            alpha = log(alpha./(1-alpha));
            gammaP = [1;1] ;
%             gammaP = [1;3] ;
            gammaP = log(gammaP);
            gammaR = [1;1] ;
            
%             gammaR = [1;3] ;
            gammaR = log(gammaR);

            for iP=1:2
                theta{iP} = [alpha(iP,:) , gammaP(iP,:) , gammaR(iP,:)];
            end
            priors.muTheta = theta{iP};
            priors.SigmaTheta = 0*eye(dim.n_theta);

            % Initial conditions
                    x0 = zeros(dim.n,1);
                    priors.muX0 = x0;
                    priors.SigmaX0 = zeros(dim.n);
                % Hyperparameters
                    % State noise precision (only for dynamical systems) 
                    alpha = Inf;
                    % Measurement noise precision 
                    sigma = Inf;

% 2/ algorithmic  options
options.priors          = priors;   % include priors in options structure
options.inG             = inG;      % input structure (grid)
options.inF             = inF;      % input structure (grid)
options.DisplayWin      = 0;
options.GnFigs          = 1;        % disable annoying figures
options.dim             = dim;
options.verbose         = 1;
options.checkGrads      = 0;
options.TolFun          = 2e-4;
options.binomial        = 1;
options.multisession.split = [];
options.multisession.split = repmat( trialNumber ,1, blockNumber); 
options.multisession.fixed.theta = 1:dim.n_theta ; % same parameters for every blocks
options.multisession.fixed.phi = 1:dim.n_phi ;
% options.multisession.fixed.X0 = 1:dim.n ;
options.skipf = zeros(1,dim.n_t);
% options.skipf(initialTrialLoss)= 1;options.skipf(initialTrialGain)= 1; % identity transition [ x1=x0 ] , instead of evolution function [ x1 = f(x0)] for initial conditions;

    

% 3/ model simulation 
    % simulation loops
    nIteration = 40; % number of iteration (usefull for monte carlo simulations)
    for iD = [1] % design loop
        for iM = [1] % model loop 
            for iP = [1:2]  % parameter loop
                for iteration = 1:nIteration
             

                    % experimental factor
    block = NaN(1,trialNumber*blockNumber);
    for iBlock = 1:blockNumber
        block(1 + trialNumber*(iBlock-1):trialNumber*(iBlock))= iBlock*ones(1,trialNumber);
    end
    trial = repmat([1:trialNumber],1,blockNumber);
    pairValence  = repmat([randperm(trialNumber/2) randperm(trialNumber/2)],1,blockNumber); % valence state    
    lottery=ones(1,n_t);
    for i=[1:3,(trialNumber/4+1:trialNumber/4+3)]
        lottery(pairValence==i)=-1;                                             % -1 = unlikely outcome (25%), % 1 = likely outcome (75%)
    end
    pairValence = -2*mod(pairValence-1,2)+1;
            
    % initial conditions
    initialTrialLoss =[];initialTrialGain=[];finalTrialLoss=[];finalTrialGain=[];
    for iBlock = 1:blockNumber
        initialTrialLoss(iBlock) = find(pairValence(block==iBlock) == -1,1,'first') + (iBlock-1)*trialNumber;
        initialTrialGain(iBlock) = find(pairValence(block==iBlock) == 1,1,'first') + (iBlock-1)*trialNumber;
        finalTrialLoss(iBlock) = find(pairValence(block==iBlock) == -1,1,'last') + (iBlock-1)*trialNumber;
        finalTrialGain(iBlock) = find(pairValence(block==iBlock) == 1,1,'last') + (iBlock-1)*trialNumber;
    end

    previousValence = zeros(1,trialNumber*blockNumber);
    previousValence(setdiff(find(trial),min(initialTrialLoss,initialTrialGain)))  = pairValence( ~ismember(find(trial),max(finalTrialLoss,finalTrialGain)));
    previousChoice = zeros(1,trialNumber*blockNumber);
    previousOutcome = zeros(1,trialNumber*blockNumber);

    bin=zeros(size(trial));
    nTrialByValence=numel(trial)/(numel(unique(block))*2);
    nBin=6;
    binSize = round(nTrialByValence/nBin);
    nBin=nTrialByValence/binSize;
    bin(pairValence==-1)=repmat(sort(repmat(1:nBin, 1, binSize)),1, numel(unique(block)));
    bin(pairValence==1)=repmat(sort(repmat(1:nBin, 1, binSize)),1, numel(unique(block)));

    u = [  pairValence ; previousValence ; previousChoice ; previousOutcome ];    
    
    fb.inH.lottery = lottery;
    fb.inH.pairValence = pairValence;
    fb.h_fname = h_fname;        
    fb.indy = 3;
    fb.indfb = 4;
    
    options.skipf(initialTrialLoss)= 1;options.skipf(initialTrialGain)= 1; % identity transition [ x1=x0 ] , instead of evolution function [ x1 = f(x0)] for initial conditions;

                    
                    
                    % predictions
                    [y,x,x0,eta,e,u] = simulateNLSS_fb(n_t,f_fname,g_fname,theta{iP},phi{iP},u,alpha,sigma,options,x0,fb);

%                     % inversion
%                         % prior densities
%                             % model parameters
%                                 % Observation parameters 
%                                 priors.muPhi = zeros(dim.n_phi,1);
%                                 priors.SigmaPhi = 0*eye(dim.n_phi);
%                                 [ priors.muPhi(1) , priors.SigmaPhi(1,1) ] = norm2log_norm(3,3,'reverse'); % comment: eg. Kr ( constraint: [0 Inf[ with a log-normal prior density == exponential transform )
%                                 % Evolution parameters 
%                                 priors.muTheta = zeros(dim.n_theta,1);
%                                 priors.SigmaTheta = 1.7*eye(dim.n_theta); % comment: eg. learning rate
%                                 % Initial conditions
%                                 priors.muX0 = zeros(dim.n,1);
%                                 priors.SigmaX0 = zeros(dim.n);
%                             % hyperparameters
%                                 % State noise precision (only for dynamical systems)  
%                                 [priors.a_alpha,priors.b_alpha]=getHyperpriors(nanvar(Y(1,:)),0.6,0.99) ; % proportion of expected variance predicted by the model
%                                 % Measurement noise precision 
%                                 sigma = Inf;
% 
%                     [posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
%                     out.diagnostics     = VBA_getDiagnostics(posterior,out);
%                         % parametric estimates
%                         inversion.parameters.posterior = posterior;
%                         inversion.parameters.Kr = exp(posterior.muPhi(1)); % inverse transform for the mode of posterior density
% 
%                     % statistical efficiency
%                     [efficiency] = VBA_designEfficiency(f_fname,g_fname,dim,options,u,'parameters'); % statistical power in parametric estimation
% %                     [efficiency] = VBA_designEfficiency(f_fname,g_fname,dim,options,u,'models'); % statistical power in model comparison
                    
                    % extract
%                     simulation.predictions.observations(iD,iM,iP,iteration) = y;
%                     routine_descriptiveStat(y)
%                     
%                     simulation.identifiability.posteriorCovariance(iD,iM,iP,iteration) = out.diagnostics.C; % == confusion matrix ;
%                     simulation.identifiability.estimationError(iD,iM,iP,iteration).phi(1) =  exp(posterior.muPhi(1)) - phi(1) ; % == estimation error ;
%                     simulation.statEfficiency.identifiability(iD,iM,iP,iteration) = efficiency;
% %                     simulation.statEfficiency.divergence(iD,iM,iP,iteration) = efficiency;
%                     simulation.design.instance(iD,iM,iP) = u;
%                     simulation.design.notation(iD,iM,iP) = designList{iD};
%                     simulation.model.instance(iD,iM,iP).specifications = specifications;
%                     simulation.model.instance(iD,iM,iP).g_function = g_fname{iM};
%                     simulation.model.instance(iD,iM,iP).f_function = f_fname{iM};
%                     simulation.model.notation(iD,iM,iP) = modelList{iM};
                    simulation(iM,iP).parameters.instance.phi = phi{iP};
                    simulation(iM,iP).parameters.instance.theta = theta{iP};
%                     simulation.parameters.notation(iD,iM,iP) = paramList{iP};
                    
                     simulation(iM,iP).output.optimalChoice(iteration,:) = y;
            simulation(iM,iP).output.probOptimalChoice(iteration,:) = y-e ;
            simulation(iM,iP).output.qValues(iteration,:,:) = x ;
            simulation(iM,iP).input.pairValence(iteration,:) = u(1,:) ;
            simulation(iM,iP).input.previousChoice(iteration,:) = u(3,:) ;
            simulation(iM,iP).input.previousOutcome(iteration,:) = u(4,:);
            simulation(iM,iP).input.outcome(iteration,:) = [u(4,2:end) 0 ];
            simulation(iM,iP).input.block(iteration,:)  = block;
            simulation(iM,iP).input.trial(iteration,:)  = trial;
            simulation(iM,iP).input.bin(iteration,:)  = bin;
            learningCurve = nanmean(tools.tapply(y-e, {u(1,:), bin, block}, @mean),3);
            simulation(iM,iP).output.negativeLearningCurve(iteration,:)  = learningCurve(1,:);
            simulation(iM,iP).output.positiveLearningCurve(iteration,:)  = learningCurve(2,:); 
            
            repetitionChoice = NaN(1,numel(y)); previousOutcome = NaN(1,numel(y));
            valCounter = 0; valNames = {'Loss';'Gain'};
            for iVal = unique(pairValence)
                valCounter = valCounter + 1;
             eval([ 'repetitionChoice(u(1,:)==iVal & ~ismember(find(trial),initialTrial' valNames{valCounter} ')) = (y(u(1,:)==iVal & ~ismember(find(trial),initialTrial' valNames{valCounter} ')) == y(u(1,:)==iVal & ~ismember(find(trial),finalTrial' valNames{valCounter} ')));' ]);
             eval([ 'previousOutcome(u(1,:)  == iVal & ~ismember(find(trial),initialTrial' valNames{valCounter} ')) =  simulation(iM,iP).input.outcome(iteration,u(1,:) ==iVal & ~ismember(find(trial),finalTrial' valNames{valCounter} '));' ]);
            end

            simulation(iM,iP).output.repetitionChoiceByOutcome(iteration,:)  = nanmean(tools.tapply(repetitionChoice, {ones(1,numel(repetitionChoice)),previousOutcome, block}, @mean),3);


                    
                    
                end
            end
        end
    end
    
% 4/ display simulation 
% 
%  hf = figure('color',[1 1 1]);
%             ha = axes('parent',hf,'nextplot','add');
%             plot(ha,mean(simulation(iM).output.optimalChoice,1),'kx')
%             plot(ha,mean(simulation(iM).output.probOptimalChoice,1),'r')
%             legend(ha,{'y: agent''s choices','p(y=1|theta,phi,m): behavioural tendency'})

            % Fig1 : Learning Curves %
            hf = figure('color',[1 1 1]);
            hf.Units = 'centimeters'; hf.Position = [ 1 1 35 15];
            
            ha = axes('parent',hf,'nextplot','add');
            ls = {'-','--'};
            paramValue = {'alpha = 0.10','alpha = 0.50'};
            hold on;
            subplot(1,2,1);
            hold on;
            for iP = 1:2
                [h(iP)] = plot(1:numel(learningCurve(1,:)),[mean(simulation(iM,iP).output.positiveLearningCurve,1)],'g');
%                 [h(iP),hp] = boundedline(1:numel(learningCurve(1,:)),[mean(simulation(iM,iP).output.positiveLearningCurve,1)],[sem(simulation(iM,iP).output.positiveLearningCurve,1)],'g', 'alpha');
                set(h(iP),'Color','g','LineWidth',2,'LineStyle',ls{iP});%set(hp,'FaceColor','g');
                h3=scatter(1:numel(learningCurve(1,:)), mean(simulation(iM,iP).output.positiveLearningCurve,1) ,80,'g','filled');
            end
            
            for proba = [0.5 0.66]
                plot([0:numel(learningCurve(1,:))+1],ones(numel(learningCurve(1,:))+2,1)*proba,'Color',[0.6 0.6 0.6],'LineWidth',1,'LineStyle','--');
            end
            
            % axis properties
            set(gca,'XTick',[1:numel(learningCurve(1,:))]);
            xlabel('blocs d''essais','FontSize',18);
            ylabel('Proportion de choix corrects','FontSize',18);
            title('Apprentissage aux gains','FontSize',18);
            hleg=legend([h(1) h(2)],texlabel(paramValue{1}),texlabel(paramValue{2}),'location','SouthWest');set(hleg,'FontSize',14);

            ylim([0 1]);xlim([1 6]);
            
            
            subplot(1,2,2);
            hold on;
            for iP = 1:2
                [h(iP)] = plot(1:numel(learningCurve(1,:)),[mean(simulation(iM,iP).output.negativeLearningCurve,1)],'r');
%                 [h(iP),hp] = boundedline(1:numel(learningCurve(1,:)),[mean(simulation(iM,iP).output.negativeLearningCurve,1)],[sem(simulation(iM,iP).output.negativeLearningCurve,1)],'r', 'alpha');
                set(h(iP),'Color','r','LineWidth',2,'LineStyle',ls{iP});%set(hp,'FaceColor','r');
                h3=scatter(1:numel(learningCurve(1,:)), mean(simulation(iM,iP).output.negativeLearningCurve,1) ,80,'r','filled');
            end
            for proba = [0.5 0.83]
                plot([0:numel(learningCurve(1,:))+1],ones(numel(learningCurve(1,:))+2,1)*proba,'Color',[0.6 0.6 0.6],'LineWidth',1,'LineStyle','--');
            end
            % axis properties
            set(gca,'XTick',[1:numel(learningCurve(1,:))]);
            xlabel('blocs d''essais','FontSize',18);
            ylabel('Proportion de choix corrects','FontSize',18);
            title('Apprentissage aux pertes','FontSize',18);
            hleg=legend([h(1) h(2)],texlabel(paramValue{1}),texlabel(paramValue{2}),'location','SouthWest');set(hleg,'FontSize',14);

            ylim([0 1]);xlim([1 6]);

%             % Fig2 : Repetition Plot %
%             hf = figure('color',[1 1 1]);
%             ha = axes('parent',hf,'nextplot','add');
%             hold on;
%             
%             hBar = bar([ mean(simulation(iM).output.repetitionChoiceByOutcome(:,1),1) ,...
%                          mean(simulation(iM).output.repetitionChoiceByOutcome(:,2),1) ,...
%                          mean(simulation(iM).output.repetitionChoiceByOutcome(:,3),1) ]);
%             hBar.FaceColor =  [0.7 0.7 0.7];
%             
%             h = errorbar([1:3],...
%                         [ mean(simulation(iM).output.repetitionChoiceByOutcome(:,1),1) ,...
%                          mean(simulation(iM).output.repetitionChoiceByOutcome(:,2),1) ,...
%                          mean(simulation(iM).output.repetitionChoiceByOutcome(:,3),1) ],...
%                    [sem(simulation(iM).output.repetitionChoiceByOutcome(:,1),1),...
%                     sem(simulation(iM).output.repetitionChoiceByOutcome(:,2),1),...
%                     sem(simulation(iM).output.repetitionChoiceByOutcome(:,3),1)]);
%                 h.LineStyle = 'none'; 
%     
%             for proba = [0.5]
%                 plot([0:4],ones(5,1)*proba,'Color',[0.6 0.6 0.6],'LineWidth',1,'LineStyle','--');
%             end
%             % axis properties
%             set(gca,'XTick',[1:3]);
%             set(gca,'XTickLabel',{'perte','neutre','gain'},'FontSize',18);
%             xlabel('conséquence précédente','FontSize',18);
%             ylabel('Proportion de choix répétés','FontSize',18);
%             ylim([0 1]);

    
