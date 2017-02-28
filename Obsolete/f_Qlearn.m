function  [fx] = f_Qlearn(x,P,u,inF)
% evolution function: quality-value delta rule
% [fx] = f_Qlearn(x,P,u,inF)
% This function solve a (2 states, 2 actions) reinforcement learning problem
%
% INPUT
%   - x_t : Q-values (for states = {loss,gain} and actions = {0,1}) (4x1)
%   - P : learning-rates (inverse logistic) ((1 + inF.stateDependantLR + inF.outcomeWeight)x1)
%     * inF.stateDependantLR==0
%       * P: learning-rate
%     * inF.stateDependantLR==1
%       * P(1): loss-dependant learning-rates
%       * P(2): gain-dependant learning-rates
%     * inF.outcomeWeight==1
%       * P: augmented with P(end+1) = outcome weight
%     * inF.outcomeWeight==2
%       * P: augmented with P(end+1) = loss weight & P(end+2) = gain weight
%   - u : 
%     * u(1) = valence(t) {loss,gain}={-1,1} (transformed into states = {1,2} )
%     * u(2) = lastValence = valence(t-1)
%     * u(3) = lastAction  = action(t-1) {0,1} 
%     * u(4) = lastOutcome = outcome(t-1) {-1,0,1} 
%   - inF : 
%     * counterfactual: option to specify the counterfactual constraint on Q-values {0,1} 
%     * stateDependantLR: option to specify the state-dependancy on learning rate {0,1} 
%     * outcomeWeight: option to specify the weight of outcome in the updating equation {0,1,2} 
%     * extraFactor = @(alpha,@q_update): homemade function handle to adapt the @f_Qlearn function to an extra-factor
% OUTPUT
%   - fx: updated Q-values (for states = {loss,gain} and actions = {0,1}) (4x1)

%default argument
% if ~exist('inF.extraFactor')
%     inF.extraFactor = @(alpha,q_update) [alpha,q_update];
% end

%%----------inputs--------------%%

valence = u(1);
state = (valence+3)/2; % [-1,1] ==> [1,2]

lastValence = u(2);
lastState = (lastValence+3)/2; % [-1,1] ==> [1,2]
lastAction = u(3)+1; % [0,1] ==> [1,2]
counterfactualLastAction = (1-u(3))+1; % [0,1] ==> [2,1]
lastOutcome = u(4);
%  lastFeedback = u(4)+2; % [-1,0,1] ==> [1,2,3]


%%----------parameters--------------%%

switch inF.outcomeWeight
    case 0
        alpha = 1./(1+exp(-P)); % learning rates are bounded between 0 and 1.
        gamma = [1,1,1]; % outcome weight fixed at 1
    case 1
        alpha = 1./(1+exp(-P(end-1))); % learning rates are bounded between 0 and 1.
        gamma = repmat(exp(P(end)),1,3); % outcome weight constrained in [0 : +Inf ]
    case 2
        alpha = 1./(1+exp(-P(end-2))); % learning rates are bounded between 0 and 1.
        gamma = [ exp(P(end-1)), 1, exp(P(end)) ]; % outcome weight constrained in [0 : +Inf ] 
end

switch inF.stateDependantLR
    case 0
        alpha=repmat(alpha,1,2);
end



%%------------evolution function: quality-value delta rule--------------%%

q(1,:)=x(1:2); q(2,:)=x(3:4);
fq=q;

% updating action value
update  = @(Q_t,s_t,a_t,fb_t,gamma,alpha) Q_t(s_t,a_t) + alpha(s_t)*(gamma(fb_t+2)*fb_t - Q_t(s_t,a_t));

fq(lastState,lastAction) = update(q,lastState,lastAction,lastOutcome,gamma,alpha);


% updating counterfactual action value   
if inF.counterfactual==1   
           fq(lastState,counterfactualLastAction) = - fq(lastState,lastAction);  
end


fx(1:2)=fq(1,:);  fx(3:4)=fq(2,:);


end
