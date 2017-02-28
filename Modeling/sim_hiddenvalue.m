function [ prediction ] = sim_hiddenvalue( param )


%% initialize
mu = param.mu;
sigma = param.sigma;
alpha = param.alpha;
mu2 = 0;
sigma2=1;

% functions
    % rank-expected hidden value
        % parameters
%         a=(mu^2)/sigma;
%         b = (mu)/sigma;
%         logmu = log(mu);
        % distribution
    %     gdist = makedist('Gamma','a',a,'b',b);
%         gdist = makedist('Lognormal','mu',logmu,'sigma',sigma);
        gdist = makedist('normal','mu',mu,'sigma',sigma);
        gdist2 = makedist('normal','mu',mu2,'sigma',sigma2);

    nq = 24; 
    ev = nan(nq,1);
    ev2 = nan(nq,1);
    for iq= 1:nq
        interval = icdf(gdist,[(iq-1)/nq iq/nq]);
        tdist = truncate(gdist,interval(1),interval(2));
        ev(iq) = mean(tdist);
        
        interval = icdf(gdist2,[(iq-1)/nq iq/nq]);
        tdist = truncate(gdist2,interval(1),interval(2));
        ev2(iq) = mean(tdist);
    end
    
    % sigmoid transform
    sigchoice = @(x) 1./(1+exp(-x));
    sigrating = @(x) 2.*sigchoice(alpha.*x)-1;


% inputs
    % rating
    seq = randperm(nq);

    % choice-1D
    choice_pair = [];
    seq1 = [ [nq/2:-1:1] ; [nq/2+1:1:nq] ];
    seq2 = [ [nq/4:-1:1],[3*nq/4:-1:nq/2+1] ; [nq/4+1:1:nq/2],[3*nq/4+1:1:nq] ];
    seq3 = [ [1:2:nq-1] ; [2:2:nq] ];
    seq4 = [ [1:1:nq/2] ; [nq/2+1:1:nq] ];
    choice_pair = [seq1,seq2,seq3,seq4];
    side = [ones(1,nq),2*ones(1,nq)];
    ind = randperm(nq*2);
    side = side(ind);
    choice_left = choice_pair(side+[0:2:(nq*4-1)]);
    choice_right = choice_pair((3-side)+[0:2:(nq*4-1)]);
    choice_pair = [choice_left ; choice_right];
    
    % choice-2D
    choice_pair = [];
    seq1 =  [ [1:nq] ; [1:nq]  ];
    seq2 =  [ [1:nq] ; [nq:-1:1]  ];
    choice_pair = [seq1,seq2];

    
%% predictions
% rating
rating = sigrating(ev(seq));

% choice-1D
dv = ev(choice_right)-ev(choice_left);
rchoice = sigchoice(dv);
rpredicted = mean(rchoice(dv>=0));
lpredicted = mean(1-rchoice(dv<0));
predictedChoice = mean([lpredicted,rpredicted]);

% choice-2D
dv = ev(choice_pair(1,:))+ev2(choice_pair(2,:));
gochoice = sigchoice(dv);
acceptanceRate = mean(gochoice);

%% output
prediction = struct;
prediction.rating = rating;
prediction.predictedChoice = predictedChoice;
prediction.acceptanceRate = acceptanceRate;

end

