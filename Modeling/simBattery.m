function [ y , U ] = simBattery( param )



            % inputs 
            % items
            task = [ 1.*ones(1,24*2) , 2.*ones(1,48*2) , 3.*ones(1,48*1) , 4.*ones(1,48*2) ] ;
            
            dimension =  [[ 1.*ones(1,24) , 4.*ones(1,24) ],...
                          [ 1.*ones(1,48) , 4.*ones(1,48) ],...
                          [ 8.*ones(1,48) ],...
                          [ 1.*ones(1,48) , 4.*ones(1,48) ]];

            x1 = [[ randperm(24) , randperm(24) ],nan(1,48*5)];
            
            x2 = nan(1,48*5+24*2); x3 = nan(1,48*5+24*2);
            [~,ind1] = sort(param.V1);             [~,ind2] = sort(param.V2);
            x2(24*2+1:24*2+48) = [ ind1(13:24) , ind1(6:-1:1) , ind1(18:-1:13) , ind1(2:13) , ind1(1:12) ];
            x3(24*2+1:24*2+48) = [ ind1(12:-1:1) , ind1(7:12) , ind1(19:24) , ind1(1:12) , ind1(13:24)];
            x2(96+1:96+48) = [ ind2(13:24) , ind2(6:-1:1) , ind2(18:-1:13) , ind2(2:13) , ind2(1:12) ];
            x3(96+1:96+48) = [ ind2(12:-1:1) , ind2(7:12) , ind2(19:24) , ind2(1:12) , ind2(13:24)];            

            x4 = nan(1,48*5+24*2); x5 = nan(1,48*5+24*2);
            x4(144+1:144+48) =  [ ind1(1:24) , ind1(1:24) ];
            x5(144+1:144+48) =  [ ind2(1:24) , ind2(24:-1:1) ];
            
            x6 = nan(1,48*5+24*2); x7 = nan(1,48*5+24*2);
            x6(192+1:192+48) =  [ ind1(randperm(12)+12) , ind1(randperm(12)+12) , ind1(randperm(12)+12) ,  ind1(randperm(12)+12) ];
            x7(192+1:192+48) =  [ ind1(randperm(12)) ,  ind1(randperm(12)) , ind1(randperm(12)) , ind1(randperm(12)) ];
            x6(240+1:240+48) =  [ ind2(randperm(12)+12) , ind2(randperm(12)+12) , ind2(randperm(12)+12) ,  ind2(randperm(12)+12) ];
            x7(240+1:240+48) =  [ ind2(randperm(12)) ,  ind2(randperm(12)) , ind2(randperm(12)) , ind2(randperm(12)) ];
            
            w1 = zeros(1,48*5+24*2);
            w2 = zeros(1,48*5+24*2);
            w3 = zeros(1,48*5+24*2);
            w4 = zeros(1,48*5+24*2);
            w5 = zeros(1,48*5+24*2);
            w6 = zeros(1,48*5+24*2);
            w7 = zeros(1,48*5+24*2);
            delays=([1 3 7 14 21 30 45 100 365 1000])./1000;
            
            w8 = nan(1,48*5+24*2);
            w8(192+1:192+96) = [ randi([1,10],1,96) ];
            w8(192+1:192+96) = delays(w8(192+1:192+96));

            s3 = [randi([0,1],1,288)*2-1] ;
            s4 = [randi([0,1],1,288)*2-1] ;
            
            
            U = [task ; dimension ;
                x1 ; x2 ; x3 ; x4 ; x5 ; x6 ; x7 ;
                w1 ; w2 ; w3 ; w4 ; w5 ; w6 ; w7 ; w8 ;
                s3 ; s4 ];
            
            U(isnan(U)) = 0;

            
        % model structure 
        f_name  = @f_evolution_battery;
        g_fname = @g_observation_battery; % observation function (delay discounted value)

        % model dimension
        dim = struct('n',4,... % number of hidden states
            'n_u',numel(U(1,:)),...  % number of time samples or trials
            'n_theta',1 ,...   % number of evolution parameters
            'n_phi',15 + 24*4 + 15 ,...     % number of observation parameters
            'p',4,...          % output (data) dimension
            'n_t',numel(U(1,:)));   % number of time samples or trials
        
        %% ITERATION ACROSS MODELS
        % initialize

                    
            phi = [ ones(1,15),...
                   (0.5).*ones(1,24*4),...
                   [1 0 0 0 0   1   0 0 0   0 0 0 0   0 0 ]];
            phi([12,15]) = [param.kRD , param.kED];
            phi([15+1:15+24]) = param.V1;
            phi([15+24*3+1:15+24*3+24]) = param.V2;
            phi([15+24*4+1]) = param.alpha;
            phi([15+24*4+2]) = param.bR;
            phi([15+24*4+5]) = param.bE;
            phi([15+24*4+6]) = param.beta;
            phi([15+24*4+14]) = param.bm;
            theta = 0;
            
            param.type = repmat({'Phi'},1,dim.n_phi);
            param.transform.direct = [ repmat({@identity},1,15) ,...
                                       repmat({@identity},1,24*4) ,...
                                       repmat({@identity},1,6) , repmat({@identity},1,9)  ];
            inG.transform = param.transform.direct(ismember(param.type,'Phi'));

            sigma = param.eta ;                % measurement noise precision


        %-------------- options ---------------% 
        inG.modelName         = 'dimensionSample';
        options.inG             = inG;      % input structure (grid)
        options.dim=dim;
        options.extended=1;
        sources(1) = struct('out',1,'type',0); % BOLD signal (gaussian, dim=4)
        sources(2) = struct('out',2:4,'type',1); % first binary response
        options.sources=sources;

        %----------------- VBA Inversion ----------------%
        [y,x,x0,eta,e] = simulateNLSS(dim.n_t,f_name,g_fname,theta,phi,U,[],sigma,options);


end

