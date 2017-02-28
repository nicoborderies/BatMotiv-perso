%%
clc;
clear all;
% close all;
%%
alpha = 1;
mu = -0.5;
sigma = 5;
t0 = 0.4;
theta = 2;

value = exp( randn(1,72).*sigma ) + mu;
sig1 = @(x) sigpos(alpha.*x,1,'direct');
sig1 = @(x) sig(alpha.*x-1);
sig2 = @(x) sig(x);
shanon = @(p) -(p.*log(p) + (1-p).*(log(1-p))); 
uncertainty = @(dv) exp(t0 + theta.*shanon(sig2(dv)));

rating = sig1(value);
ind1 = randperm(72);
ind2 = randperm(72);
dr = value(ind2) - value(ind1);
choice = sig2(dr);
choice(dr<0) = 1 - choice(dr<0);
rt = uncertainty(dr);

bins = [-5:0.1:5];
col = [1 1 1]*0.5;

%%
figure;
histogram(value,bins,'Normalization','Probability','EdgeColor','none','FaceColor',col);
xlabel('hidden value distribution');


figure;
%     % sigpos function 
%         [x,y] = fplot(sig1,[ bins(1) bins(end) ],'-k');hold on;
%         plot(x,y,'-k');
% 
%         alpha = 0.5;
%         sig1 = @(x) sigpos(alpha.*x,1,'direct');
%         [x,y] = fplot(sig1,[ bins(1) bins(end) ],'-k');hold on;
%         plot(x,y,'-k');
%     
    % sigmoid function 
        [x,y] = fplot(sig1,[ bins(1) bins(end) ],'-k');hold on;
        plot(x,y,'-k');

        alpha =2;
        sig1 = @(x) sig(alpha.*x-1);
        [x,y] = fplot(sig1,[ bins(1) bins(end) ],'-k');hold on;
        plot(x,y,'-k');

xlabel('hidden value');
ylabel('rating value');


figure;hold on;
histogram(rating,[0:0.05:1],'Normalization','Probability','EdgeColor','none','FaceColor',col);
xlabel('rating distribution');


figure;
% sigmoid function 
[x,y] = fplot(sig2,[ -5 5],'-k');hold on;
plot(x,y,'-k');
xlabel('V2 - V1');
ylabel('P(choice=2)');

figure;
% uncertainty function 
[x,y] = fplot(uncertainty,[ -10 10],'-k');hold on;
plot(x,y,'-k');
xlabel('V2 - V1');
ylabel('uncertainty(choice)');

figure;hold on;
b = barplot(1,mean(choice),std(choice)*0,col);
ax = gca;
% ax.XTick = [1];
% ax.XTickLabel = {'choice(1D) accuracy'};
ax.YLim = [0 1];
xlabel('choice(1D) accuracy');

figure;hold on;
histogram(rt,[0:0.1:5],'Normalization','Probability','EdgeColor','none','FaceColor',col);
xlabel('response time distribution');
