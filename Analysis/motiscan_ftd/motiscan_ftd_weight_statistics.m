%% motiscan_ftd_weight_statistics

%%  statistics
%%% 1.  acceptance rate comparisons
%%%% rewards
clc;fprintf('\n');
y1 = stat.weight.acceptRate(stat.rating.group==1);
y2 = stat.weight.acceptRate(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median acceptance rate, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

%%% 2. regression coeficients comparisons
%%%% rewards
clc;fprintf('\n');
y1 = stat.weight.kR(stat.rating.group==1);
y2 = stat.weight.kR(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median reward weight, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

%%%% efforts
clc;fprintf('\n');
y1 = stat.weight.kE(stat.rating.group==1);
y2 = stat.weight.kE(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median effort weight, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

    
%%% 3. comparison of coefficients to null value 
%%%% rewards
clc;fprintf('\n');
y2 = stat.weight.kR(stat.rating.group==2);
%%% sign test
[p,h,stats] = signtest(y2,0);
n = numel(y2);
npos = stats.sign;
fprintf('sign test against 0, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);
%%%%  efforts
y2 = stat.weight.kE(stat.rating.group==2);
%%% sign test
[p,h,stats] = signtest(y2,0);
n = numel(y2);
npos = stats.sign;
fprintf('sign test against 0, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);

%%% 4. comparison of dimension aggregation 
%%%% proportion of joint significance
clc;fprintf('\n');
y1 = stat.weight.significant_RE(stat.rating.group==1);
y2 = stat.weight.significant_RE(stat.rating.group==2);
%%%% medians
mu = nanmean([y1,y2]);
fprintf('mean proportion of significant joint weights , control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% chi-2 test
[h,p,chi,df] = chi2test2(y1,y2);
fprintf('Chi-2 test, ? = %.3g, p = %.2e \n',chi,p);

%%%% minimal weight
fprintf('\n');
y1 = stat.weight.min_abs_kR_kE(stat.rating.group==1);
y2 = stat.weight.min_abs_kR_kE(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median minimal absolute weight, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%% 5. comparison of dimension aggregation to chance level
%%%% minimal weight
fprintf('\n');
y1 = stat.weight.min_signed_kR_kE(stat.rating.group==1);
y2 = stat.weight.min_signed_kR_kE(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median minimal signed weight, bvFTD = %.2g \n',mu(2));
%%% sign test
[p,h,stats] = signtest(y2,0);
n = numel(y2);
npos = stats.sign;
fprintf('sign test against 0, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);

%%% 6. comparison of relative reward/effort weight ratio
fprintf('\n');
y1 = stat.weight.kR(stat.rating.group==1) + stat.weight.kE(stat.rating.group==1);
y2 = stat.weight.kR(stat.rating.group==2) + stat.weight.kE(stat.rating.group==2);
% y1 = stat.weight.kR(stat.rating.group==1)./ ([ stat.weight.kR(stat.rating.group==1) - stat.weight.kE(stat.rating.group==1) ]);
% y2 = stat.weight.kR(stat.rating.group==2)./ ([ stat.weight.kR(stat.rating.group==2) - stat.weight.kE(stat.rating.group==2) ]);

%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median sum of reward/effort weights, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);


