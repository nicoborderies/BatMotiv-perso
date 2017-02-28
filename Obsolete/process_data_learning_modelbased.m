function [result] = process_data_learning_modelbased(dataStructure, option)
% process_data_learning_modelbased process the observations collected with
% the learning task of the MBB Battery, processing is a model-based
% analysis
%
% [result] = process_data_learning_modelbased(dataStructure, option)
%
% INPUT
%       - data structure from the load_data_learning function
%       - option is an optionnal structure with the following fields
%           * sessionList : session to process (default option is to process all sessions) (1xnSession)
%           * metadata : information about the subjects (eg. treatment, diagnostic, demographic informations)
%              * iSub : subject identification number
%           * model : structure containing specifications concerning the models implemented, ordered by their identification number
%              * counterfactual : option to specify the counterfactual constraint on Q-values {0,1} 
%              * stateDependantLR : option to specify the state-dependency on learning rate {0,1} 
%              * stateDependantT : option to specify the state-dependancy on temperature {0,1} 
%           * extraFactor = @(option.metadata,U,extraFactorF,extraFactorG,dim) : function handle that adapt the mapping of the model according to extra factors contained in the metadata structure (eg. a group factor)
%
% OUTPUT
%       - result containing the following field
%       -------------------------------------------------------------------------------------------------------------
%           * inversion : structure containing models implemented, ordered by their identification number
%               * posterior: posterior estimates for every model's parameter
%               * out: informations about model's inversion
%               * fit: fitting metrics about the model
%               * parameters: model's parameters
%               * specifications: specifications concerning the model implemented
%
%
% NB : sessionNumber is converted to get the real number of SESSIONS (i.e. days)
%    : blockNumber is the number of block by session

%default arguments
if nargin<2;
    option.nBin=6;
    option.sessionList=unique(dataStructure.condition.sessionNumber); 
elseif ~isfield(option,'sessionList')
    option.sessionList=unique(dataStructure.condition.sessionNumber); 
elseif ~isfield(option,'extraFactor')
    if~isfield(option,'metadata');metadata=[];end;
    extraFactor = @(metadata,U,f_fname,g_fname,dim) [U,f_fname,g_fname,dim];
end


% Filter the processing with respect to the sessionList
blockNumber     = dataStructure.condition.blockNumber(ismember(dataStructure.condition.sessionNumber,option.sessionList));
isOptimalChoice = dataStructure.behavior.isOptimalChoice(ismember(dataStructure.condition.sessionNumber,option.sessionList));
pairValence     = dataStructure.condition.pairValence(ismember(dataStructure.condition.sessionNumber,option.sessionList));
outcome         = dataStructure.behavior.outcome(ismember(dataStructure.condition.sessionNumber,option.sessionList));



% Get the real "session" number
sessionNumber = floor((dataStructure.condition.sessionNumber-1)/max(dataStructure.condition.blockNumber))+1;
sessionNumber = sessionNumber(ismember(dataStructure.condition.sessionNumber,option.sessionList));
nSession=max(sessionNumber);
nBlock=max(dataStructure.condition.blockNumber);

% Metadata
% if option.metadata.TREATMENT(1,option.metadata.groupid(option.metadata.iSub))=='P'
%     treatment   =sessionNumber;
% elseif option.metadata.TREATMENT(1,option.metadata.groupid(option.metadata.iSub))=='C'
%     treatment   =[sessionNumber==2; sessionNumber==1];
% end

%Rename and subset to non-neutral pair
sessionNumber   = sessionNumber(pairValence~=0);
% treatment       = treatment (dataStructure.condition.pairValence~=0);
blockNumber     = blockNumber(pairValence~=0);
isOptimalChoice = isOptimalChoice(pairValence~=0);
pairValence     = pairValence(pairValence~=0);
outcome = outcome(pairValence~=0);
blockSize = numel(isOptimalChoice(blockNumber(sessionNumber==1)==1));


%% ANALYSES
% ===========================================================================

for iModel = 1:numel(option.model) % Loop across models
    for iSession = unique(sessionNumber) % Loop across session
        
        fprintf('processing learning_modelbased: subject (%d/%d), model (%d/%d), session(%d/%d) \n',...
            option.metadata.iSub,numel(option.metadata.groupid),iModel,numel(option.model),iSession,numel(unique(sessionNumber)));

        % ------------- Mappings ----------------- %    
        % outputs 
        Y = [isOptimalChoice(sessionNumber==iSession)];                    

        % inputs 
        inG.pairValence = pairValence(sessionNumber==iSession);

            % initial conditions
            nNeg=[];nPos=[];jBlock=0;
            for iBlock = unique(blockNumber(sessionNumber==iSession))
                jBlock = jBlock+1;
                nNeg(iBlock) = find(inG.pairValence(blockNumber(sessionNumber==iSession)==iBlock) == -1,1,'first') + (jBlock-1)*blockSize;
                nPos(iBlock) = find(inG.pairValence(blockNumber(sessionNumber==iSession)==iBlock) == 1,1,'first') + (jBlock-1)*blockSize;
            end

        lastPairValence = pairValence(sessionNumber==iSession); lastPairValence = [0 lastPairValence(1:end-1)];
        lastOptimalChoice = isOptimalChoice(sessionNumber==iSession) ; lastOptimalChoice = [0 lastOptimalChoice(1:end-1)];
        lastOutcome = outcome(sessionNumber==iSession);lastOutcome = [0 lastOutcome(1:end-1)];

        U = [  pairValence(sessionNumber==iSession) ; lastPairValence ; lastOptimalChoice ; lastOutcome];

        % model structure 
        f_fname = @f_Qlearn; % evolution function (Q-value learning)
        g_fname = @g_softmax; % observation function (softmax selection)

        % options
        inF.counterfactual = option.model{iModel}.counterfactual;
        inF.stateDependantLR = option.model{iModel}.stateDependantLR;
        inG.stateDependantT = option.model{iModel}.stateDependantT;


        % model space
        dim = struct('n',4,... % number of hidden states
            'n_theta',1+inF.stateDependantLR ,...    % number of evolution parameters
            'n_phi',1+inG.stateDependantT,...      % number of observation parameters
            'p',1,...          % output (data) dimension
            'n_t',numel(Y));   % number of time samples or trials
        
        % extra-factor mapping
%         [U,inF.extraFactor,inG.extraFactor,dim] = extraFactor(option.metadata,U,extraFactorF,extraFactorG,dim);


        %-------------- priors definition ---------------% 
        % Model parameters
            % Observation parameters 
            priors.muPhi = zeros(dim.n_phi,1);
            priors.SigmaPhi = 1*eye(dim.n_phi);
            % Evolution parameters 
            priors.muTheta = zeros(dim.n_theta,1);
            priors.SigmaTheta = 1.7*eye(dim.n_theta);
            % Initial conditions
            priors.muX0 = zeros(dim.n,1);
            priors.SigmaX0 = zeros(dim.n);
        % Hyperparameters
            % State noise precision (only for dynamical systems) 
            priors.a_alpha = Inf;
            priors.b_alpha = 0;


        %-------------- options ---------------% 

        % inversion scheme
        options.priors          = priors;   % include priors in options structure
        options.inG             = inG;      % input structure (grid)
        options.inF             = inF;      % input structure (grid)
        options.DisplayWin      = 1;
        options.GnFigs          = 0;        % disable annoying figures
        options.dim             = dim;
        options.verbose         = 1;
        options.checkGrads      = 0;
        options.TolFun          = 2e-4;
        %options.isYout          =  zeros(1,numel(Y));   
        options.binomial        = 1;
        options.multisession.split = [];
        for iBlock=unique(blockNumber(sessionNumber==iSession))
           options.multisession.split = [options.multisession.split  blockSize]; 
        end
        options.multisession.fixed.theta = 1:dim.n_theta ; % same parameters for every blocks
        options.multisession.fixed.phi = 1:dim.n_phi ;
        options.skipf = zeros(1,dim.n_t);
        for iBlock=unique(blockNumber(sessionNumber==iSession))
            options.skipf(nNeg(iBlock))= 1;options.skipf(nPos(iBlock))= 1; % identity transition [ x1=x0 ] , instead of evolution function [ x1 = f(x0)] for initial conditions;
        end

        %----------------- VBA Inversion ----------------%
        [posterior,out] = VBA_NLStateSpaceModel(Y,U,f_fname,g_fname,dim,options);
        out.diagnostics     = VBA_getDiagnostics(posterior,out);

        % fitting metrics
        result.inversion{iModel}.posterior{iSession} = posterior;
        result.inversion{iModel}.out{iSession} = out;
        result.inversion{iModel}.fit(iSession).R2 = out.fit.R2;
        result.inversion{iModel}.fit(iSession).BCA = out.fit.acc;
        result.inversion{iModel}.fit(iSession).bayesian_p = 1./(1+exp(out.F - out.diagnostics.LLH0));
        result.inversion{iModel}.fit(iSession).logE = out.F;
        result.inversion{iModel}.parameters{iSession}.alpha = 1./(1+exp(-posterior.muTheta));
        result.inversion{iModel}.parameters{iSession}.beta = exp(posterior.muPhi);
    end
    
    result.inversion{iModel}.specifications = option.model{iModel};
end



end