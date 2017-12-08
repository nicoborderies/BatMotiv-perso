%% motiscan_ftd_choice_statistics

%%  statistics
%%% 1. choice accuracy comparisons
%%%% rewards
clc;fprintf('\n');
y1 = stat.choice.choice_accuracy_r(stat.rating.group==1);
y2 = stat.choice.choice_accuracy_r(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median choice accuracy, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% deviance to normality
[h,p] = kstest(y1);
fprintf('kolmogorov test, p = %.2e \n',p);
[h,p] = kstest(y2);
fprintf('kolmogorov test, p = %.2e \n',p);
%%%% ttest
[h,p,~,test] = ttest2(y1,y2);
t = test.tstat;
df = test.df;
p = p;
fprintf('two-tailed unpaired t-test, t(%d) = %.3g, p = %.2e \n',df,t,p);
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%%  efforts
fprintf('\n');
y1 = stat.choice.choice_accuracy_e(stat.rating.group==1);
y2 = stat.choice.choice_accuracy_e(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median choice accuracy, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% deviance to normality
[h,p] = kstest(y1);
fprintf('kolmogorov test, p = %.2e \n',p);
[h,p] = kstest(y2);
fprintf('kolmogorov test, p = %.2e \n',p);
%%%% ttest
[h,p,~,test] = ttest2(y1,y2);
t = test.tstat;
df = test.df;
p = p;
fprintf('two-tailed unpaired t-test, t(%d) = %.3g, p = %.2e \n',df,t,p);
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%% 2. comparison to chance level
%%%% rewards
clc;fprintf('\n');
y1 = stat.choice.choice_accuracy_r(stat.rating.group==1);
y2 = stat.choice.choice_accuracy_r(stat.rating.group==2);
%%% sign test
[p,h,stats] = signtest(y2,0.5);
n = numel(y2);
npos = stats.sign;
fprintf('sign test against 0.5, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);
%%%%  efforts
y1 = stat.choice.choice_accuracy_e(stat.rating.group==1);
y2 = stat.choice.choice_accuracy_e(stat.rating.group==2);
%%% sign test
[p,h,stats] = signtest(y2,0.5);
n = numel(y2);
npos = stats.sign;
fprintf('sign test against 0.5, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);


%%% 3. logistic regression weight comparisons
%%%% rewards
clc;fprintf('\n');
y1 = stat.choice.dv_weight_r(stat.rating.group==1);
y2 = stat.choice.dv_weight_r(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median weight, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% deviance to normality
[h,p] = kstest(y1);
fprintf('kolmogorov test, p = %.2e \n',p);
[h,p] = kstest(y2);
fprintf('kolmogorov test, p = %.2e \n',p);
%%%% ttest
[h,p,~,test] = ttest2(y1,y2);
t = test.tstat;
df = test.df;
p = p;
fprintf('two-tailed unpaired t-test, t(%d) = %.3g, p = %.2e \n',df,t,p);
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%%  efforts
fprintf('\n');
y1 = stat.choice.dv_weight_e(stat.rating.group==1);
y2 = stat.choice.dv_weight_e(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median weight, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% deviance to normality
[h,p] = kstest(y1);
fprintf('kolmogorov test, p = %.2e \n',p);
[h,p] = kstest(y2);
fprintf('kolmogorov test, p = %.2e \n',p);
%%%% ttest
[h,p,~,test] = ttest2(y1,y2);
t = test.tstat;
df = test.df;
p = p;
fprintf('two-tailed unpaired t-test, t(%d) = %.3g, p = %.2e \n',df,t,p);
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);
    