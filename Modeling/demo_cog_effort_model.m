% demo_cog_effort_model

% parameters
beta = [2 1];
kG = 1;
kC = 3;

% inputs
cong = -1;
x = [0 1];
y = [0 1];
gain = 1;

% equations
    % default response proba
    q = sig( beta(1)*cong*abs(x(2)-x(1)) + beta(2)*abs(y(2)-y(1)) );
    % controlled response proba
    p = @(u) sig( (1/(1+u))*beta(1)*cong*abs(x(2)-x(1)) + (1+u)*beta(2)*abs(y(2)-y(1)) ) ; 
    % information cost
    cost = @(u) KLDiv( [p(u),1-p(u)] , [q,1-q] ) ; 
    % adaptation cost
%     cost = @(u) sqrt( (beta(1)/(1+u) - beta(1))^2 + (beta(2)*(1+u) - beta(2))^2 );
    % value function
    value = @(u) kG*gain*p(u) - kC*cost(u);
    negvalue = @(u) -value(u);

    
% predictions
[u_s,v_s] = fminbnd(negvalue,0,10);
v_s = -v_s;
pcor = p(u_s);

% plots
    % value function
    figure;
    hold on;
    fplot(p,[0 10]);
    fplot(cost,[0 10]);
    fplot(value,[0 10]);