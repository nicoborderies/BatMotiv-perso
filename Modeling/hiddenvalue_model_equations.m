%% hiddenvalue_model_equations
%

%% reset
clc;
clear all;
% close all;


%% parameters

alpha = 1;
mu = 0;
sigma = 1;

bm = 0;
bp = 0;
kd = 1;

t0 = 0.2;
theta = 1;


%% functions

valdensity = @(x) (1/sqrt(2*pi*sigma^2)).*exp(-(x-mu).^2/(2*sigma^2) );
valquantile = @(x) mu + sigma*sqrt(2)*erfinv(2*x-1);

sigchoice = @(x) 1./(1+exp(-x));
sigrating = @(x) 2.*sigchoice(alpha.*x)-1;
% probarating = @(x) 

discount = @(x,d) x./(1+kd*d);

shanon = @(p) -(p.*log(p) + (1-p).*(log(1-p))); 
uncertainty = @(dv) exp(t0 + theta.*shanon(sigchoice(dv)));



%% displays

% 1-hidden value density 
f1 = figure; hold on;
for mu = [ 0 2]
    for sigma = [1 2]
        valdensity = @(x) (1/sqrt(2*pi*sigma^2)).*exp(-(x-mu).^2/(2*sigma^2) );
        fplot(valdensity,[-5 5]);
    end
end
xlabel('hidden value');
ylabel('probability density');
title('hidden value density');
box on;
axis([-5 5 0 0.5]);
legend({'\mu = 0, \sigma = 1',...
       '\mu = 0, \sigma = 2',...
       '\mu = 2, \sigma = 1',...
       '\mu = 2, \sigma = 2'});
legend('boxoff');
fig.Position = [ 100 100 500 400];
set_all_properties('FontName','Arial Narrow','FontWeight','normal','FontSize',16,...
                    'LineWidth',1.5,'Interpreter','tex');

%% 2-quantile-expected value 
f2 = figure; hold on;
for mu = [ 0 2]
    for sigma = [1 2]
        valquantile = @(x) mu + sigma*sqrt(2)*erfinv(2*x-1);
        fplot(valquantile,[0 1]);
    end
end
xlabel('quantile rank');
ylabel('expected value');
title('quantile-expected value');
box on;
axis([0 1 -5 5]);
legend({'\mu = 0, \sigma = 1',...
       '\mu = 0, \sigma = 2',...
       '\mu = 2, \sigma = 1',...
       '\mu = 2, \sigma = 2'});
legend('boxoff');
fig.Position = [ 100 100 500 400];
set_all_properties('FontName','Arial Narrow','FontWeight','normal','FontSize',16,...
                    'LineWidth',1.5,'Interpreter','tex');
                
%% 3-rating prediction
f3 = figure; hold on;
for alpha = [ 1 2]
    sigrating = @(x) 2.*sigchoice(alpha.*x)-1;
    fplot(sigrating);
end
xlabel('hidden value');
ylabel('rating value');
title('rating transform');
box on;
axis([-5 5 -1 1]);
legend({'\alpha = 1',...
       '\alpha = 2'});
legend('boxoff');
fig.Position = [ 100 100 500 400];
set_all_properties('FontName','Arial Narrow','FontWeight','normal','FontSize',16,...
                    'LineWidth',1.5,'Interpreter','tex');
                
%% 4-choice prediction
f4 = figure; hold on;
fplot(sigchoice,'k');
xlabel('\Delta value');
ylabel('P( choice = max(V) )');
title('choice transform');
box on;
axis([-5 5 0 1]);
legend('boxoff');
fig.Position = [ 100 100 500 400];
set_all_properties('FontName','Arial Narrow','FontWeight','normal','FontSize',16,...
                    'LineWidth',1.5,'Interpreter','tex');
                
%% 5-response time prediction
f5 = figure; hold on;
fplot(uncertainty,'k');
xlabel('\Delta value');
ylabel('response time');
title('uncertainty transform');
box on;
axis([-5 5 0 3]);
legend('boxoff');
fig.Position = [ 100 100 500 400];
set_all_properties('FontName','Arial Narrow','FontWeight','normal','FontSize',16,...
                    'LineWidth',1.5,'Interpreter','tex');
