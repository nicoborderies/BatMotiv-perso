%% gompertz_function_simulations

% close all;
% clc;
% clear;

% parameters
a = 1;
b = 1;
c = 1;
d = 0;

% function
gomp = @(x)  a.*exp(-b.*exp(-c.*x))+d ;


% inputs
xx = rand(1,1000);


% display 
f1 = figure; hold on;
fplot(gomp,[-5 5]);

% hist
f2 = figure; hold on;
histogram(gomp(xx),20);

