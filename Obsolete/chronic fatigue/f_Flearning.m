function  [fx] = f_Flearning(x,P,u,inF)
% evolution function: quality-value delta rule
% [fx] = f_Flearning(x,P,u,inF)
% This function solve a 3 hidden force levels learning problem
%
% INPUT
%   - x_t : F-values (3x1)
%   - P : learning-rates (inverse logistic) (1x1)
%       * P(1): learning-rate
%   - u : 
%     * u(1) = level(t) {30% , 50% , 70%} (transformed into states = {1,2,3} )
%     * u(2) = previousLevel {30% , 50% , 70%} (transformed into states = {1,2,3} )= level(t-1)
%     * u(3) = previousOutcome {too little , correct , too much} = outcome(t-1) {1,0,-1} 
%   - inF : 
%
%
%
% OUTPUT
%   - fx: updated F-values (for levels = {30% , 50% , 70%}) (3x1)


%%----------inputs--------------%%

level = u(1);
 % [30,50,70] ==> [1,2,3]

previousLevel = u(2);
% [30,50,70] ==> [1,2,3]

previousOutcome = u(3);




%%----------parameters--------------%%

alpha = 1./(1+exp(-P)); % learning rates are bounded between 0 and 1.




%%------------evolution function: quality-value delta rule--------------%%

F_t_1 = x;

F_t = F_t_1;

% updating action value
update  = @(F_t,l_t,o_t,alpha) F_t(l_t) + alpha*(o_t);

F_t(level) = update(F_t_1,previousLevel,previousOutcome,alpha);


fx=F_t;

end