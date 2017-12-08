%% motiscan_ftd_grip_statistics

%%  statistics
%%% 1.  mean fpeak comparisons
clc;fprintf('\n');
y1 = stat.grip.mean_fpeak(stat.rating.group==1);
y2 = stat.grip.mean_fpeak(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median force peak average, control = %.4g N, bvFTD = %.4g N \n',mu(1),mu(2));
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

clc;fprintf('\n');
y1 = stat.grip.mean_nfpeak(stat.rating.group==1);
y2 = stat.grip.mean_nfpeak(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median force peak intercept, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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


%%% 2.  incentive coeficient comparisons
clc;fprintf('\n');
y1 = stat.grip.incentive_fpeak(stat.rating.group==1);
y2 = stat.grip.incentive_fpeak(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median incentive coeficient, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

%%% 3.  ntrial coeficient comparisons
clc;fprintf('\n');
y1 = stat.grip.ntrial_fpeak(stat.rating.group==1);
y2 = stat.grip.ntrial_fpeak(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median trial number coeficient, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

%%% 4. maximal force peak comparisons
clc;fprintf('\n');
y1 = stat.grip.max_fpeak(stat.rating.group==1);
y2 = stat.grip.max_fpeak(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median maximal force peak, control = %.4g N, bvFTD = %.4g N \n',mu(1),mu(2));
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

%%% 5. calibration force peak comparisons
clc;fprintf('\n');
y1 = stat.grip.calib_fpeak(stat.rating.group==1);
y2 = stat.grip.calib_fpeak(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median maximal force peak, control = %.4g N, bvFTD = %.4g N \n',mu(1),mu(2));
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

%%% 6.  response time comparisons
clc;fprintf('\n');
y1 = stat.grip.rt(stat.rating.group==1);
y2 = stat.grip.rt(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median response time, control = %.3g s, bvFTD = %.3g s \n',mu(1),mu(2));
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


%% 
%%% 7. parameters comparisons
%%%% kr
clc;fprintf('\n');
y1 = stat.grip.kr(stat.rating.group==1);
y2 = stat.grip.kr(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median kr, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% ke
fprintf('\n');
y1 = stat.grip.ke(stat.rating.group==1);
y2 = stat.grip.ke(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median ke, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% kr
fprintf('\n');
y1 = stat.grip.kf(stat.rating.group==1);
y2 = stat.grip.kf(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median kf, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% tau
fprintf('\n');
y1 = stat.grip.tau(stat.rating.group==1);
y2 = stat.grip.tau(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median tau, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% R2_force
fprintf('\n');
y1 = stat.grip.R2_force(stat.rating.group==1);
y2 = stat.grip.R2_force(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median R2_force, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% R2_velocity
fprintf('\n');
y1 = stat.grip.R2_velocity(stat.rating.group==1);
y2 = stat.grip.R2_velocity(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median R2_velocity, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);


%% 
%%% 8. oxymeter parameters comparisons
%%%% mean hr response
clc;fprintf('\n');
y1 = stat.grip.mean_hr_response(stat.rating.group==1);
y2 = stat.grip.mean_hr_response(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median HR response, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% mean incentive hr response
fprintf('\n');
y1 = stat.grip.incentive_hr_response(stat.rating.group==1);
y2 = stat.grip.incentive_hr_response(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median incentive HR response, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% mean ntrial hr response
fprintf('\n');
y1 = stat.grip.ntrial_hr_response(stat.rating.group==1);
y2 = stat.grip.ntrial_hr_response(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median ntrial HR response, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% mean corrected incentive hr response
fprintf('\n');
y1 = stat.grip.incentive_hr_response_corrected(stat.rating.group==1);
y2 = stat.grip.incentive_hr_response_corrected(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median corrected incentive HR response, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% mean corrected ntrial hr response
fprintf('\n');
y1 = stat.grip.ntrial_hr_response_corrected(stat.rating.group==1);
y2 = stat.grip.ntrial_hr_response_corrected(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median corrected ntrial HR response, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%% 9. oxymetric-behavioral correlation of incentive markers
fprintf('\n');
% within-group (ftd)
% y1 = stat.grip.incentive_hr_response(stat.rating.group==2);
% y2 = stat.grip.incentive_fpeak(stat.rating.group==2);
% across-group
y1 = stat.grip.incentive_hr_response;
y2 = stat.grip.incentive_fpeak;
%%%% spearman correlation
[rho,p] = corr(y1,y2,'type','Spearman','rows','pairwise');
fprintf('Spearman correlation two-tailed test, rho = %.3g, p = %.2e \n',rho,p);

