%% motiscan_ftd_learning_statistics

%%  statistics
%%% 1. correct rate comparisons
clc;fprintf('\n');
y1 = nanmean([stat.learning.correct_gain(stat.rating.group==1),...
              stat.learning.correct_loss(stat.rating.group==1)  ],2);
y2 = nanmean([stat.learning.correct_gain(stat.rating.group==2),...
              stat.learning.correct_loss(stat.rating.group==2)  ],2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median correct rate, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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


%%% 1. correct rate comparison to chance level
clc;fprintf('\n');
%%% gain
y1 = stat.learning.correct_gain(stat.rating.group==1);
y2 = stat.learning.correct_gain(stat.rating.group==2);
%%% sign test
[p,h,stats] = signtest(y2,0.5);
n = numel(y2);
npos = stats.sign;
fprintf('sign test against 0.5, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);
%%%%  efforts
fprintf('\n');
y1 = stat.learning.correct_loss(stat.rating.group==1);
y2 = stat.learning.correct_loss(stat.rating.group==2);
%%% sign test
[p,h,stats] = signtest(y2,0.5);
n = numel(y2);
npos = stats.sign;
fprintf('sign test against 0.5, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);

% 2. valence decomposition
clc;fprintf('\n');
y = [stat.learning.correct_loss;...
      stat.learning.correct_gain ];
x = nominal([ ones(nsub,1) ; 2*ones(nsub,1)]) ;
g = nominal([ stat.rating.group ; stat.rating.group]) ; 

m = fitlm(table(x,g,y),'y ~ 1 + x*g ');
t = anova(m);df1 = t.DF(1);df2 = t.DF(end);
mu = splitapply(@nanmedian,y,double(x));
disp(t);
fprintf('ANOVA, group, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(1),t.pValue(1));
fprintf('ANOVA, dimension, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(2),t.pValue(2));
fprintf('ANOVA, group*dimension, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(3),t.pValue(3));
fprintf('median loss correct rate = %.2g \n',mu(1));
fprintf('median gain correct rate = %.2g \n',mu(2));

%%% 3. information weight comparisons
%%%% gains
clc;fprintf('\n');
y1 = stat.learning.lor_weight_gain(stat.rating.group==1);
y2 = stat.learning.lor_weight_gain(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median reliable information weight|gain, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

%%%% loss
fprintf('\n');
y1 = stat.learning.lor_weight_loss(stat.rating.group==1);
y2 = stat.learning.lor_weight_loss(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median reliable information weight|loss, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

%%% 3. information weight comparison to chance level
clc;fprintf('\n');
%%% gain
y1 = stat.learning.lor_weight_gain(stat.rating.group==1);
y2 = stat.learning.lor_weight_gain(stat.rating.group==2);
%%% sign test
[p,h,stats] = signtest(y2,0);
n = numel(y2);
npos = stats.sign;
fprintf('sign test against 0, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);
%%%%  efforts
fprintf('\n');
y1 = stat.learning.lor_weight_loss(stat.rating.group==1);
y2 = stat.learning.lor_weight_loss(stat.rating.group==2);
%%% sign test
[p,h,stats] = signtest(y2,0);
n = numel(y2);
npos = stats.sign;
fprintf('sign test against 0, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);



%%% 4. delay weight comparisons
%%%% gains
clc;fprintf('\n');
y1 = stat.learning.delay_weigth(stat.rating.group==1);
y2 = stat.learning.delay_weigth(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median delay weight, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

%%% 5. repetitive choice comparison
%%%% neutral
clc;fprintf('\n');
y1 = stat.learning.mean_repchoice(stat.rating.group==1);
y2 = stat.learning.mean_repchoice(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median repetitive choice rate, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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
%%% 6. parameters comparisons
%%%% alpha
clc;fprintf('\n');
y1 = stat.learning.alpha(stat.rating.group==1);
y2 = stat.learning.alpha(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median alpha, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% kr
fprintf('\n');
y1 = stat.learning.kr(stat.rating.group==1);
y2 = stat.learning.kr(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median kr, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% kp
fprintf('\n');
y1 = stat.learning.kp(stat.rating.group==1);
y2 = stat.learning.kp(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median kp, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% bm
fprintf('\n');
y1 = stat.learning.bm(stat.rating.group==1);
y2 = stat.learning.bm(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median bm, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% bp
fprintf('\n');
y1 = stat.learning.bp(stat.rating.group==1);
y2 = stat.learning.bp(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median bp, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

%%%% BCA
fprintf('\n');
y1 = stat.learning.BCA(stat.rating.group==1);
y2 = stat.learning.BCA(stat.rating.group==2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median classification accuracy, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);