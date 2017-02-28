%% hiddenvalue_simulations
%
% nicolas borderies - 12/2016

clc;
clear all;
close all;

%% initialize

% model parameters
% murange = [0.1:0.5:3.1];
% sigmarange = [0.1:0.5:3.1];

% murange = [1];
% sigmarange = [0.25,1,4];

% murange = [0,1,2];
% sigmarange = [1];

% murange2 = [0];
% sigmarange2 = [1];

murange = [1];
sigmarange = [1,1,1];
alpharange = [1,2,4];

% paramvalues = combnk(prange,2);
paramvalues = combvec(murange,sigmarange)';
MU = paramvalues(:,1);
SIGMA = paramvalues(:,2);
ALPHA = alpharange;

% sim parameters
nsim = numel(MU);

% predictions
nq = 24;
predictedChoice = nan(1,nsim);
predictedRating = nan(nq,nsim);
acceptanceRate = nan(1,nsim);

%% simulations


for isim = 1:nsim
    
    fprintf(' simulation: %d/%d \n',isim,nsim);

    param.mu = MU(isim);
    param.sigma =  SIGMA(isim);
    param.alpha =  ALPHA(isim);

    [ prediction ] = sim_hiddenvalue( param );
    
    predictedChoice(isim) = prediction.predictedChoice;
    predictedRating(:,isim) = prediction.rating;
    acceptanceRate(isim) = prediction.acceptanceRate;

end

%% display

% x = [0:0.1:10];
% gpdf = pdf(gdist,x);
% % figure;
% hold on;
% plot(x,gpdf);

% % figure;
% hold on;
% plot([1:nq],ev);

% 

figure; hold on;
for isim = 1:nsim
     x = [0:0.05:1];
%     [f,b] = histcounts(predictedRating(:,isim),'Normalization','probability','BinWidth',0.05);
    [f,b] = ksdensity(predictedRating(:,isim),[0:0.01:1]);
    plot(b,f);
end
ylabel('probability density (%)');
title('rating distribution');
% h = legend({'mu=0','mu=1','mu=2'});h.Box='off';
% h = legend({'sigma=0.25','sigma=1','sigma=4'});h.Box='off';
h = legend({'alpha=1','alpha=2','alpha=4'});h.Box='off';


% figure;
% [X,Y] = meshgrid(murange,sigmarange);
% X(:) = MU;
% Y(:) = SIGMA;
% Z = nan(size(X));
% Z(:) = predictedChoice;
% surf(X,Y,Z);
% xlabel('mu');
% ylabel('sigma');
% zlabel('predicted choice (%)');
% zlim([0.5 1]);

% figure;
% [X,Y] = meshgrid(murange,sigmarange);
% X(:) = MU;
% Y(:) = SIGMA;
% Z = nan(size(X));
% Z(:) = acceptanceRate;
% surf(X,Y,Z);
% xlabel('mu');
% ylabel('sigma');
% zlabel('acceptance rate (%)');
% zlim([0.5 1]);

