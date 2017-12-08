%% motiscan_ftd_clinical_stat

% 
subset1 = (subtab.group=='CONTROL');
subset2 = (subtab.group=='FTD');

% sex ratio
%%% proportion of females
clc;fprintf('\n');
y1 = double(subtab.sex(subset1)=='F');
y2 = double(subtab.sex(subset2)=='F');
%%%% medians
mu = sum([y1,y2]);
n = size([y1,y2],1);
fprintf('female/male ratio , control = %.2g/%.2g, bvFTD = = %.2g/%.2g \n',mu(1),n-mu(1),mu(2),n-mu(2));
%%%% chi-2 test
[h,p,chi,df] = chi2test2(y1,y2);
fprintf('Chi-2 test, ?(1) = %.3g, p = %.2e \n',chi,p);

% age
clc;fprintf('\n');
y1 = subtab.age(subset1);
y2 = subtab.age(subset2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median age, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

% education level
clc;fprintf('\n');
y1 = subtab.educationLevel(subset1);
y2 = subtab.educationLevel(subset2);
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median education level, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

% apathy score
clc;fprintf('\n');
y1 = subtab.aes(subset2);
y2 = subtab.hes(subset2);
%%%% medians
mu = sum([y1,y2]>=14);
fprintf('N(Starkstein AES>=14) = %.3g \n',mu(1));
fprintf('N(Starkstein HES>=14) = %.3g \n',mu(2));

% anosognosia
%%%% medians
mu = nanmedian([y1,y2]);
fprintf('median apathy score, auto = %.2g, hetero = %.2g \n',mu(1),mu(2));
%%%% wilcoxon rank-sum test
[p,h,test] = ranksum(y1,y2);
U = test.ranksum;
fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);

% cognitive decline score
clc;fprintf('\n');
y2 = subtab.Mattis(subset2);
%%%% medians
mu = nanmedian([y2]);
fprintf('median Mattis score = %.3g \n',mu);
%%% sign test
[p,h,stats] = signtest(y2,123,'tail','left');
n = numel(y2);
npos = stats.sign;
fprintf('left-tailed sign test against 123, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);

% dysexecutive score
fprintf('\n');
y2 = subtab.FAB(subset2);
%%%% medians
mu = nanmedian([y2]);
fprintf('median FAB score = %.3g \n',mu);
%%% sign test
[p,h,stats] = signtest(y2,12,'tail','left');
n = numel(y2);
npos = stats.sign;
fprintf('left-tailed sign test against 12, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);

fprintf('\n');
y2 = subtab.HAD(subset2);
%%%% medians
mu = nanmedian([y2]);
fprintf('median HAD score = %.3g \n',mu);
%%% sign test
[p,h,stats] = signtest(y2,13,'tail','right');
n = numel(y2);
npos = stats.sign;
fprintf('right-tailed sign test against 13, npos = %.3g, n = %.3g, p = %.2e \n',npos,n,p);

%%%% display table
mean_stat = varfun(@nanmean,subtab(:,[1 8 10 11 12 13 14 15 16 17 18 19]),'GroupingVariables','group');
sd_stat = varfun(@nanstd,subtab(:,[1 8 10 11 12 13 14 15 16 17 18 19]),'GroupingVariables','group');

groupstat = [ mean_stat , sd_stat(:,[3:end]) ];
groupstat = groupstat(:,[1 2 3 3+11 4 4+11 5 5+11 6 6+11 7 7+11 8 8+11 9 9+11 10 10+11 11 11+11 12 12+11 13 13+11]);
count_female = @(x) sum(x=='F');
sex_stat = varfun(count_female,subtab(:,[1 9]),'GroupingVariables','group');

disp(groupstat);
