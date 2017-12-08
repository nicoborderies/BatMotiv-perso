%% motiscan_ftd_rating_statistics

%% 3.3 statistics
%%% 3.3.1 average rating comparisons
%%%% rewards
fprintf('\n');
y1 = stat.rating.mean_r(stat.rating.group==1);
y2 = stat.rating.mean_r(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median rating average, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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
y1 = -stat.rating.mean_e(stat.rating.group==1);
y2 = -stat.rating.mean_e(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median rating average, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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

%%% 3.3.2 paired rating comparisons
%%%% rewards
fprintf('\n');
y1 = median(itemRatings(stat.rating.group==1,:,1));
y2 = median(itemRatings(stat.rating.group==2,:,1));
y = double(y2>y1);
mu = mean(y);
[h,p,ci,z] = ztest(y,0.5,std(y));
fprintf('fraction of items = %.2g\n',mu);
fprintf('z-test, z = %.3g, p = %.2e \n',z,p);
%%%% efforts
fprintf('\n');
y1 = median(itemRatings(stat.rating.group==1,:,4));
y2 = median(itemRatings(stat.rating.group==2,:,4));
y = double(y2>y1);
mu = mean(y);
[h,p,ci,z] = ztest(y,0.5,std(y));
fprintf('fraction of items = %.2g\n',mu);
fprintf('z-test, z = %.3g, p = %.2e \n',z,p);

%%% 3.3.3 subdomain comparisons
%%%% rewards
fprintf('\n');
y = [ stat.rating.mean_food ; stat.rating.mean_good ] ;
x = nominal([ ones(nsub,1) ; 2*ones(nsub,1)]) ;
g = nominal([ stat.rating.group ; stat.rating.group]) ; 
x = ([ ones(nsub,1) ; 2*ones(nsub,1)]) ;

[p,tbl,stats,terms] = anovan(y,{x g},'model','interaction','varnames',{'x','g'},'display','off');
m = fitlm(table(x,g,y),'y ~ -1 + x*g ');
t = anova(m);df1 = t.DF(1);df2 = t.DF(end);
fprintf('ANOVA, group*subdomain, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(3),t.pValue(3));
c = multcompare(stats,'Dimension',[1 2],'CType','lsd');
fprintf('Least Significant Difference procedure, bvFTD-control|food_rating, p = %.2e \n',c(2,end));
fprintf('Least Significant Difference procedure, bvFTD-control|good_rating, p = %.2e \n',c(5,end));
m = fitlm(table(x,g,y),'y ~ 1 + x*g ');
t = anova(m);df1 = t.DF(1);df2 = t.DF(end);
fprintf('ANOVA, group*food_rating, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(2),t.pValue(2));
fprintf('ANOVA, group*good_rating, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(3),t.pValue(3));

%%%% efforts
clc;fprintf('\n');
y = [ stat.rating.mean_motor ; stat.rating.mean_cognitive ] ;
x = nominal([ ones(nsub,1) ; 2*ones(nsub,1)]) ;
g = nominal([ stat.rating.group ; stat.rating.group]) ; 
x = ([ ones(nsub,1) ; 2*ones(nsub,1)]) ;

[p,tbl,stats,terms] = anovan(y,{x g},'model','interaction','varnames',{'x','g'},'display','off');
m = fitlm(table(x,g,y),'y ~ -1 + x*g ');
t = anova(m);df1 = t.DF(1);df2 = t.DF(end);
fprintf('ANOVA, group*subdomain, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(3),t.pValue(3));
c = multcompare(stats,'Dimension',[1 2],'CType','lsd');
fprintf('Least Significant Difference procedure, bvFTD-control|food_rating, p = %.2e \n',c(2,end));
fprintf('Least Significant Difference procedure, bvFTD-control|good_rating, p = %.2e \n',c(5,end));
m = fitlm(table(x,g,y),'y ~ 1 + x*g ');
t = anova(m);df1 = t.DF(1);df2 = t.DF(end);
fprintf('ANOVA, group*food_rating, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(2),t.pValue(2));
fprintf('ANOVA, group*good_rating, F(%d,%d) = %.3g, p = %.2e \n',df1,df2,t.F(3),t.pValue(3));

%%% 3.3.4 extreme rating comparisons
%%%% rewards
fprintf('\n');
y1 = stat.rating.extrem_r(stat.rating.group==1);
y2 = stat.rating.extrem_r(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median extreme rating proportion, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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
y1 = stat.rating.extrem_e(stat.rating.group==1);
y2 = stat.rating.extrem_e(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median extreme rating proportion, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
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


%%% differentiating the 2 boundaries
clc;fprintf('\n');
%%%% rewards 
y1 = itemRatings(stat.rating.group==1,:,1);
y2 = itemRatings(stat.rating.group==2,:,1);
minY1 = min(y1,[],2);maxY1 = max(y1,[],2);
minY2 = min(y2,[],2);maxY2 = max(y2,[],2);
% minimal rating
    f1 =  mean((y1==minY1),2);
    f2 =  mean((y2==minY2),2);
    mu = median([f1,f2]);
    fprintf('median minimal rating proportion, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
    %%%% wilcoxon rank-sum test
    [p,h,test] = ranksum(f1,f2);
    U = test.ranksum;
    fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);
% maximal rating
    f1 =  mean((y1==maxY1),2);
    f2 =  mean((y2==maxY2),2);
    mu = median([f1,f2]);
    fprintf('median maximal rating proportion, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
    %%%% wilcoxon rank-sum test
    [p,h,test] = ranksum(f1,f2);
    U = test.ranksum;
    fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);
    
clc;fprintf('\n');
%%%% efforts 
y1 = itemRatings(stat.rating.group==1,:,4);
y2 = itemRatings(stat.rating.group==2,:,4);
minY1 = min(y1,[],2);maxY1 = max(y1,[],2);
minY2 = min(y2,[],2);maxY2 = max(y2,[],2);
% minimal rating
    f1 =  mean((y1==minY1),2);
    f2 =  mean((y2==minY2),2);
    mu = median([f1,f2]);
    fprintf('median minimal rating proportion, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
    %%%% wilcoxon rank-sum test
    [p,h,test] = ranksum(f1,f2);
    U = test.ranksum;
    fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);
% maximal rating
    f1 =  mean((y1==maxY1),2);
    f2 =  mean((y2==maxY2),2);
    mu = median([f1,f2]);
    fprintf('median maximal rating proportion, control = %.2g, bvFTD = %.2g \n',mu(1),mu(2));
    %%%% wilcoxon rank-sum test
    [p,h,test] = ranksum(f1,f2);
    U = test.ranksum;
    fprintf('Mann-Whitney-Wilcoxon test, U = %.3g, p = %.2e \n',U,p);
    
    
%% 3.3.5 response time 
%%%% rewards
fprintf('\n');
y1 = stat.rating.responsetime_r(stat.rating.group==1);
y2 = stat.rating.responsetime_r(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median response time, control = %.4g s, bvFTD = %.4g s \n',mu(1),mu(2));
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
y1 = stat.rating.responsetime_e(stat.rating.group==1);
y2 = stat.rating.responsetime_e(stat.rating.group==2);
%%%% medians
mu = median([y1,y2]);
fprintf('median response time, control = %.4g s, bvFTD = %.4g s \n',mu(1),mu(2));
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
    