function [ result ] = predict_hiddenvalue( param )


%% inputs 
% task parameters
nq = 24;

% items
task = [ 1.*ones(1,nq) , 2.*ones(1,nq*2) , 3.*ones(1,nq*2)  ] ;

rating_rank = zeros(size(task));
rating_rank(task==1) = randperm(nq);

choice_ranks = zeros(2,numel(task));
seq1 = [ [nq/2:-1:1] ; [nq/2+1:1:nq] ];
seq2 = [ [nq/4:-1:1],[3*nq/4:-1:nq/2+1] ; [nq/4+1:1:nq/2],[3*nq/4+1:1:nq] ];
seq3 = [ [1:2:nq-1] ; [2:2:nq] ];
seq4 = [ [1:1:nq/2] ; [nq/2+1:1:nq] ];
choice_ranks(:,task==2) = [seq1,seq2,seq3,seq4];

accept_ranks = zeros(2,numel(task));
seq1 =  [ [1:nq] ; [1:nq]  ];
seq2 =  [ [1:nq] ; [nq:-1:1]  ];
accept_ranks(:,task==3) = [seq1,seq2]; % column1 = benefits, column2 = costs

inputs = [task ; rating_rank ; choice_ranks ; accept_ranks ];


%% model structure 
g_fname = @g_hiddenvalue;
options.inG.nq = nq;

% model dimension
dim = struct('n',0,... % number of hidden states
            'n_u',size(inputs,1),...  % dimension of inputs
            'n_theta',0 ,...   % number of evolution parameters
            'n_phi',4,...     % number of observation parameters
            'p',4,...          % output (data) dimension
            'n_t',numel(task));   % number of time samples or trials
        
        
% parameters
alpha = param.alpha;
mu = param.mu;
sigma = param.sigma;
theta = param.theta;
phi = [alpha , theta , mu , sigma ];

% hyperparameters
hyper_alpha = Inf;
hyper_sigma = [ param.precision , param.precision ];

% options
options.sources(1) = struct('out',1,'type',0); % normal
options.sources(2) = struct('out',2:3,'type',1); % binomial response
options.sources(3) = struct('out',4,'type',0); % normal

%%  predictions
[outputs,x,x0,eta,e] = simulateNLSS(dim.n_t,[],g_fname,[],phi,inputs,hyper_alpha,hyper_sigma,options,[]);



%% output

result = struct;
result.inputs = inputs;
result.outputs = outputs;
result.param = param;

choice = outputs(2,task==2);
% drank = double( choice_ranks(1,task==2) < choice_ranks(2,task==2) );
result.predictedChoice = mean(choice==1);
accept = outputs(3,task==3);
result.acceptRate = mean(accept==1);



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


