%% motiscan_ftd_correlations

% 
subset1 = (stat.rating.group==2);
subset2 = (subtab.group=='FTD');

% apathy score 
%%% preference parameters
fprintf('\n');
y1 = stat.battery.sd_R(subset1);
% y1 = stat.battery.sd_E(subset1);
% y1 = nanmean([stat.battery.sd_R(subset1), stat.battery.sd_E(subset1)],2);
% y1 = stat.battery.mu_R(subset1);
% y1 = stat.battery.mu_E(subset1);
y2 = subtab.aes(subset2);
% y2 = subtab.hes(subset2);
%%%% spearman correlation
[rho,p] = corr(y1,y2,'type','Spearman','rows','pairwise');
% corrplot(y1,y2);
fprintf('Spearman correlation two-tailed test, rho = %.3g, p = %.2e \n',rho,p);

%%% force parameters
fprintf('\n');
% y1 = log(stat.grip.kr(subset1));
y1 = log(stat.grip.ke(subset1));
% y1 = stat.grip.incentive_fpeak(subset1);
% y1 = stat.grip.mean_fpeak(subset1);
% y2 = subtab.aes(subset2);
y2 = subtab.hes(subset2);
%%%% spearman correlation
[rho,p] = corr(y1,y2,'type','Spearman','rows','pairwise');
% corrplot(y1,y2);
fprintf('Spearman correlation two-tailed test, rho = %.3g, p = %.2e \n',rho,p);

%%% learning parameters
fprintf('\n');
y1 = stat.learning.kr(subset1);
% y1 = stat.learning.kp(subset1);
% y1 = stat.learning.alpha(subset1);
% y1 = stat.learning.correct_gain(subset1);
% y1 = stat.learning.correct_loss(subset1);
y2 = subtab.aes(subset2);
% y2 = subtab.hes(subset2);
%%%% spearman correlation
[rho,p] = corr(y1,y2,'type','Spearman','rows','pairwise');
% corrplot(y1,y2);
fprintf('Spearman correlation two-tailed test, rho = %.3g, p = %.2e \n',rho,p);

% FBI score
fprintf('\n');
y1 = stat.battery.sd_R(subset1);
y1 = stat.battery.kD(subset1);
y1 = stat.battery.bD(subset1);
y2 = subtab.FBI(subset2);
%%%% spearman correlation
[rho,p] = corr(y1,y2,'type','Spearman','rows','pairwise');
% corrplot(y1,y2);
fprintf('Spearman correlation two-tailed test, rho = %.3g, p = %.2e \n',rho,p);

% alimentary behavior score
fprintf('\n');
y1 = stat.rating.mean_food(subset1);
% y1 = stat.rating.mean_food(subset1)-stat.rating.mean_good(subset1);
% y1 = stat.rating.extrem_r(subset1);
% y1 = stat.battery.mu_R(subset1);
% y1 = stat.battery.sd_R(subset1);
y2 = subtab.Qalim(subset2);
%%%% spearman correlation
[rho,p] = corr(y1,y2,'type','Spearman','rows','pairwise');
% corrplot(y1,y2);
fprintf('Spearman correlation two-tailed test, rho = %.3g, p = %.2e \n',rho,p);
