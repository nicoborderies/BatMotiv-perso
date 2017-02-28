function [gx] = g_Fselection(x_t,P,u_t,inG)
% observation function:    cost-benefit valuation + optimal control of
% force estimate

% [fx] = g_Fselection(x,P,u,inF)
%
% INPUT
% - x_t : Ft(levels) : force estimates at time t
% - P : Kp, Ke
% - u_t : 
%     * u(1) = level(t) {30% , 50% , 70%} (transformed into states = {1,2,3} )
% - inG :
%
%
% OUTPUT
% - gx : 
%       - gx: Scalar, predicted force


% Define inputs & parameters
%--------------------------------------------------------------------------

level           = u_t(1);



Kp          = [exp(P(1))] ;
Ke          = [exp(P(2))] ;

   
% Compute ouput
%--------------------------------------------------------------------------


gx=zeros(1,1); % intialize
    
% identity function
gx = x_t(level);


% % cost benefit valuation
% Prob = @(f) ...
% Prob = @(f) 5*exp(-(f-level).^2);
% 
% precisionCost = @(f) ...
% precisionCost = @(f) 1./(f);
% 
% motorCost = @(f) ...
% motorCost = @(f) log(1-f).^2;
% 
% value = @(f) Prob(f) - Kp*precisionCost(f) - Ke*motorCost(f);
% negValue = @(f) - value(f);
% 
% % optimisation
% gx = fminbnd(negValue,0.01,0.99);
            


end


