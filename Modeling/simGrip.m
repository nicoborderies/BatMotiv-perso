function [ y , U ] = simGrip( param )



            % inputs 
            reward=[0.01 0.2 0.5  1 5 20];
            blocks=1:10;
            trials=1:12;
            trialNumber = [1:numel(blocks)*numel(trials)];
            incentive=[];
            for nblock=blocks
                cond=[];
                for i=1:length(trials)/6
                    cond=[cond randperm(6)];
                end
                incentive=[incentive mod(cond-1,6)+1];                   % 1cents=1, 20cents=2, 50cents=3, 1euro=4, 5euros=5, 20euros=6
            end
            incentive = reward(incentive);
            valence = mod(ceil(trialNumber/12),2)*2-1 ;
            
            U = [  incentive/(max(incentive)).*100 ;
                valence ;
                trialNumber/60 ;
                nan(1,max(trialNumber)) ;
                nan(1,max(trialNumber)) ;
                repmat(param.calib,1,max(trialNumber)) ];

            
            U(isnan(U)) = 0;

            
        % model structure 
        g_fname = @g_effortSelection; % observation function (delay discounted value)

        % model dimension
        dim = struct('n',0,... % number of hidden states
            'n_u',numel(U(1,:)),...  % number of time samples or trials
            'n_theta',0 ,...   % number of evolution parameters
            'n_phi',7 ,...     % number of observation parameters
            'p',2,...          % output (data) dimension
            'n_t',numel(U(1,:)));   % number of time samples or trials
        
        %% ITERATION ACROSS MODELS
        % initialize

                    
            phi(1) = param.kR;
            phi(2) = param.kP;
            phi(3) = param.kE;
            phi(4) = param.k0;
            phi(5) = param.kF;
            phi(6) = param.tau;
            phi(7) = param.fmax - param.calib;

            
            param.type = repmat({'Phi'},1,dim.n_phi);
            param.transform.direct = [ repmat({@identity},1,7) ];
            inG.transform = param.transform.direct(ismember(param.type,'Phi'));
            sigma = [5e1 5e1] ;                % measurement noise precision


        %-------------- options ---------------% 
        inG.modelName         = 'dimensionSample';
        inG.predictYank =1;
        inG.maxObservedForce = param.calib;
        options.inG             = inG;      % input structure (grid)
        options.dim=dim;
        options.extended=1;
        sources(1) = struct('out',1,'type',0); % (gaussian)
        sources(2) = struct('out',2,'type',0); % (gaussian)
        options.sources=sources;

        %----------------- VBA Inversion ----------------%
        [y,x,x0,eta,e] = simulateNLSS(dim.n_t,[],g_fname,[],phi,U,[],sigma,options);


end

