function [estimates] = estimate_hiddenvalue(inputs,outputs,param)


%% inputs 
% task parameters
nq = 24;
task = inputs(1,:);

%% model structure 
g_fname = @g_hiddenvalue;
options.inG.nq = nq;

% model dimension
dim = struct('n',0,... % number of hidden states
            'n_u',size(inputs,1),...  % dimension of inputs
            'n_theta',0 ,...   % number of evolution parameters
            'n_phi',4,...     % number of observation parameters
            'p',4,...          % output (data) dimension
            'n_t',size(inputs,2));   % number of time samples or trials
        
        
% parameters
% phi = [alpha , theta , mu , sigma ];
priors.muPhi = ones(dim.n_phi,1); 
priors.SigmaPhi = 1e0*eye(dim.n_phi);
options.priors = priors;

% hyperparameters

% options
options.sources(1) = struct('out',1,'type',0); % normal
options.sources(2) = struct('out',2:3,'type',1); % binomial response
options.sources(3) = struct('out',4,'type',0); % normal
for i=1:2
    options.priors.a_sigma(i) = 1;
    options.priors.b_sigma(i) = 1;
end

options.isYout = zeros(size(outputs));
options.isYout(1,task~=1) = 1;
options.isYout(2,task~=2) = 1;
options.isYout(3,task~=3) = 1;
options.isYout(4,task==1) = 1;

if param.get_choice==0
    options.isYout(2,task==2) = 1;
end
if param.get_accept==0
    options.isYout(3,task==3) = 1;
end
if param.get_responsetime==0
    options.isYout(4,:) = 1;
end

options.DisplayWin = 0;


%%  inversion
[posterior,inversion] = VBA_NLStateSpaceModel(outputs,inputs,[],g_fname,dim,options);



%% output

estimates = struct;
estimates.posterior = posterior;
estimates.inversion = inversion;
estimates.param.alpha = posterior.muPhi(1);
estimates.param.theta = posterior.muPhi(2);
estimates.param.mu = posterior.muPhi(3);
estimates.param.sigma = posterior.muPhi(4);
estimates.confusionMatrix = cov2corr(posterior.SigmaPhi);


end


function [gx] = g_hiddenvalue(x_t,P,u_t,inG)

% parameters
alpha = P(1);
theta = P(2);
t_0 = 0.2;

mu = P(3);
sigma = P(4);
mu_cost = 0;
sigma_cost = 1;

% inputs
task = u_t(1) ;
rating_rank = u_t(2) ;
choice_ranks = u_t([3,4]) ;
benefit_rank = u_t(5) ;
cost_rank = u_t(6) ;

nq = inG.nq;

% functions
rank2quantile = @(x) max([ min([ (nq-x)/(nq-1) , (1-1e-3) ]) , 1e-3 ]);
valdensity = @(x) (1/sqrt(2*pi*sigma^2)).*exp(-(x-mu).^2/(2*sigma^2) );
valquantile = @(x,mu,sigma) mu + sigma*sqrt(2)*erfinv(2*x-1);
sigchoice = @(x) 1./(1+exp(-x));
sigrating = @(x) 2.*sigchoice(alpha.*x)-1;
% probarating = @(x) 
discount = @(x,d) x./(1+kd*d);
shanon = @(p) -(p.*log(p) + (1-p).*(log(1-p))); 
uncertainty = @(dv,t0,theta) exp(t0 + theta.*shanon(sigchoice(dv)));

% predictions
rating = 0;
choice = 0;
accept = 0;
responsetime = 0;

if task == 1
    ev_rating = valquantile( rank2quantile(rating_rank) ,mu,sigma);
    rating = sigrating( ev_rating );
    gx(1) = rating;
    
elseif task == 2
    dv = valquantile( rank2quantile(choice_ranks(1)),mu,sigma) - valquantile( rank2quantile(choice_ranks(2)),mu,sigma);
    choice = sigchoice( dv );
    responsetime = uncertainty(choice,t_0,theta);
    
elseif task == 3 
    ev_accept = valquantile( rank2quantile(benefit_rank),mu,sigma) + valquantile(rank2quantile(cost_rank),mu_cost,sigma_cost);
    accept = sigchoice( ev_accept );
    responsetime = uncertainty(accept,t_0,theta);

end

gx = [ rating ; choice ; accept ; responsetime ];

end


