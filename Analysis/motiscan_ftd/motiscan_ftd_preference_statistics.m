%% motiscan_ftd_preference_statistics

%%  statistics
%%% 1 parameters comparisons
%%%% mu_R
clc;fprintf('\n');
y1 = stat.battery.mu_R(stat.rating.group==1);
y2 = stat.battery.mu_R(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median mu_R, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% mu_E
fprintf('\n');
y1 = stat.battery.mu_E(stat.rating.group==1);
y2 = stat.battery.mu_E(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median mu_E, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% sigma_R
fprintf('\n');
y1 = stat.battery.sd_R(stat.rating.group==1);
y2 = stat.battery.sd_R(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median sd_R, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% sigma_E
fprintf('\n');
y1 = stat.battery.sd_E(stat.rating.group==1);
y2 = stat.battery.sd_E(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median sd_E, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% mu_R-mu_E
fprintf('\n');
y1 = stat.battery.mu_R(stat.rating.group==1) + stat.battery.mu_E(stat.rating.group==1);
y2 = stat.battery.mu_R(stat.rating.group==2) + stat.battery.mu_E(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median mu_R+mu_E, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% sigma_R/sigma_E
fprintf('\n');
y1 = stat.battery.sd_R(stat.rating.group==1)./stat.battery.sd_E(stat.rating.group==1);
y2 = stat.battery.sd_R(stat.rating.group==2)./stat.battery.sd_E(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median sd_R/sd_E, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% alpha
fprintf('\n');
y1 = stat.battery.alpha(stat.rating.group==1);
y2 = stat.battery.alpha(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median alpha, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% bD
fprintf('\n');
y1 = stat.battery.bD(stat.rating.group==1);
y2 = stat.battery.bD(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median bD, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% kD
fprintf('\n');
y1 = stat.battery.kD(stat.rating.group==1);
y2 = stat.battery.kD(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median kD, control = %.3g, bvFTD = %.3g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% fit accuracy
fprintf('\n');
y1 = stat.battery.R2_rating(stat.rating.group==1);
y2 = stat.battery.R2_rating(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median R2 rating, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

fprintf('\n');
y1 = stat.battery.BCA_choice(stat.rating.group==1);
y2 = stat.battery.BCA_choice(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median choice classification accuracy, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);


%% 2. valence decomposition of group effect on sigma
clc;fprintf('\n');
y = [stat.battery.sd_R;...
      stat.battery.sd_E ];
x = nominal([ ones(nsub,1) ; 2*ones(nsub,1)]) ;
g = nominal([ stat.rating.group ; stat.rating.group]) ; 

% robust linear model + ANOVA
m = fitlm(table(x,g,y),'y ~ 1 + x*g ','RobustOpts','on');
t = anova(m);df1 = t.DF(1);df2 = t.DF(end);
mu = splitapply(@nanmedian,y,double(x));
disp(t);
fprintf('ANOVA, group, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(1),t.pValue(1));
fprintf('ANOVA, dimension, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(2),t.pValue(2));
fprintf('ANOVA, group*dimension, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(3),t.pValue(3));
fprintf('median sd_R = %.2g \n',mu(1));
fprintf('median sd_E = %.2g \n',mu(2));
[h,p] = kstest(m.Residuals.Raw);
fprintf('ANOVA residuals; kolmogorov test, p = %.2e \n',p);

%%% 3. parameter correlations 
clc;fprintf('\n');
% within-group (ftd)
% y1 = mean( [ stat.battery.sd_R(stat.rating.group==2) , stat.battery.sd_E(stat.rating.group==2)],2);
% y2 = stat.battery.alpha(stat.rating.group==2);
% across-group
y1 = mean( [ stat.battery.sd_R , stat.battery.sd_E],2);
y2 = stat.battery.alpha;
%%%% spearman correlation
[rho,p] = corr(y1,y2,'type','Spearman','rows','pairwise');
fprintf('Spearman correlation two-tailed test, rho = %.3g, p = %.2e \n',rho,p);


