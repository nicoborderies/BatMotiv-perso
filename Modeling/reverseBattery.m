function [estimate] = reverseBattery(y, u)

     % model structure 
                f_name  = @f_evolution_battery;
                g_fname = @g_observation_battery; % observation function (delay discounted value)

            % model dimension
                 dim = struct('n',4,... % number of hidden states
                            'n_theta',1 ,...   % number of evolution parameters
                            'n_phi',15 + 24*4 + 15 ,...     % number of observation parameters
                            'p',4,...          % output (data) dimension
                            'n_t',numel(y(1,:)));   % number of time samples or trials
            
            
            % priors
                    phi = struct;
                    phi.prior.mu = [ ones(1,15),...
                                       (0).*ones(1,24*4),...
                                       [1 0 0 0 0  1   0 0 0   0 0 0 0   0 0 ],...
                                          0 ];
                    phi.prior.sigma  = [ 1.*ones(1,15),...
                                           10.*ones(1,24*4),...
                                           1.*ones(1,15),...
                                           0 ];
                    phi.prior.sigma([1:11]) = 0;
                    phi.prior.sigma([ 15 + 24*4 + 7 : 15 + 24*4 + 13 ]) = 0;
                    phi.prior.sigma([ 15 + 24*4 + 15 ]) = 0;
                    phi.prior.sigma([15+24*4+1]) = 0;

                    
                    for iDim=[2,3]
                            phi.prior.sigma(iDim)=0;
                            phi.prior.sigma(7 + iDim)=0;
                            phi.prior.sigma(11 + iDim)=0;
                            phi.prior.sigma( 15 + 24*(iDim-1) + 1 : 15 + 24*(iDim-1)  + 24)=0;
                            phi.prior.sigma(15 + 24*4 + 1 + iDim)=0;
                            phi.prior.sigma(15 + 24*4 + 13 + iDim)=0;
                    end
                    for iDim=9:10
                            phi.prior.sigma(15 + 24*4 + iDim - 1)=0;
                    end
                     phi.type = repmat({'Phi'},1,dim.n_phi);
                    phi.labels = [ {'kR2','kRI2','kP2','kE2','kR3','kP3','kE3','kR4','kRI4','kP4','kE4','kRD','kRID','kPD','kED'},...
                                    repmat({'V'},1,24*4),...
                                    {'alpha','bR','bRI','bP','bE','beta','bRE','bPE','bRP','bRD','bRID','bPD','bED','bm','bp'},...
                                    {'lambda'}];
                    item = cellfun(@num2str,num2cell([1:24]),'UniformOutput',0);
                    for idim=1:4
                        for it=1:24
                            phi.labels{15 + (idim-1)*24 + it} = [ phi.labels{15 + (idim-1)*24 + it} , num2str(idim) , '_' , num2str(it) ];
                        end
                    end
                    phi.transform.direct = [ repmat({@safepos},1,15) ,...
                                               repmat({@safepos},1,24*4) ,...
                                               { @safepos @safepos @safepos @safepos @safepos @safepos } , repmat({@identity},1,9),...
                                               {@identity} ];
                    inG.transform = phi.transform.direct(ismember(phi.type,'Phi'));

                    opt.display=0;
                   [priors]  = setParam(phi,opt);
                   
                % Initial conditions
                    priors.muX0 = ones(dim.n,1)*0;
                    priors.SigmaX0 = eye(dim.n)*0;


                    % Hyperparameters
            %             % State noise precision (only for dynamical systems) 
                            priors.a_alpha = Inf;
                            priors.b_alpha = 0;
                        % Measurement noise precision
                            [priors.a_sigma(1),priors.b_sigma(1)]=getHyperpriors(nanvar(y(1,:)),0.10,0.90) ; 
                            [priors.a_sigma(2),priors.b_sigma(2)]=getHyperpriors(nanvar(y(1,:)),0.10,0.90) ; 
      
                   
                    inG.modelName         = 'dimensionSample';
                    options.priors          = priors;   % include priors in options structure
                    options.inG             = inG;      % input structure (grid)
                    options.DisplayWin      = 1;
                    options.dim             = dim;
                    options.isYout          =  zeros(size(y,1),size(y,2));   % data exclusion
                    options.isYout(isnan(y))          =  1;   
                    for iT=1:4
                         options.isYout(iT,u(1,:)~=iT) = 1;
                    end
                    options.updateX0 = 0 ;
                    options.kernelSize  = 0 ;
                    options.checkGrads = 0 ;
                    options.extended=1;
                    options.sources(1).out  = 1;    % ratings
                    options.sources(1).type = 0;    % Normal continous data
                    options.sources(2).out  = 2:4 ;    % choices, weight, discount
                    options.sources(2).type = 1;    % binomial data

            
            [posterior,inversion] = VBA_NLStateSpaceModel(y,u,f_name,g_fname,dim,options);
            
            [phi]  = getParam(phi,posterior,opt);
            estimate = struct;
            for iP = 1:numel(posterior.muPhi)
                estimate.(phi.labels{iP}) = phi.estimate(iP);
                if phi.prior.sigma(iP)==0;
                    estimate.(phi.labels{iP}) = NaN;
                end
            end
            estimate.muV1 = mean(phi.estimate(15+1:15+24));
            estimate.sigmaV1 = std(phi.estimate(15+1:15+24));
            estimate.muV2 = mean(phi.estimate(15+72+1 : 15+72+24));
            estimate.sigmaV2 = std(phi.estimate(15+72+1 : 15+72+24));

end

