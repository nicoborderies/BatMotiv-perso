%% motiscan_opioid_crosstask_correlations
%


%% morphine effect onto force intensity & force maintenance tasks
% preparation 
y1 = nan(nsub,1);
y2 = nan(nsub,1);

% extract subject stats
for isub = 1:nsub % subject loop
    
    % stat 1
    tab = data{isub}.grip.table;
    fpeak = tab.normalizedForcePeak ;
    rt = tab.rt ;
    y = fpeak;
%     y = rt;
    trt = tab.treatment;
    trt = removecats(trt,'0');
    trt = reordercats(trt,treatmentList);
    g = double(trt);
    ycond = splitapply(@nanmean,y,g);
%     y1(isub) = ycond(3) - ycond(2);
%     y1(isub) = ycond(3) - mean(ycond);
%     y1(isub) = (ycond(3) - ycond(2)) + (ycond(1) - ycond(2)) ;
%     y1(isub) = ycond(2);
%     y1(isub) = ycond(3);
%     y1(isub) = mean(ycond);
    predictor = [ (trt=='naloxone')*(-1) + (trt=='placebo')*0 + (trt=='morphine')*1 ];
    formula = ['y ~ 1 + opioid'];
    stat = fitglm(predictor,y,formula,'VarNames',{'opioid','y'});
    coef = stat.Coefficients;
    y1(isub) = coef.Estimate(2);
    
    % stat 2
    tab = data{isub}.gripAccu.table;
    y = tab.effortRatio ;
    trt = tab.treatment;
    trt = removecats(trt,'0');
    trt = reordercats(trt,treatmentList);
    g = double(trt);
    ycond = splitapply(@nanmean,y,g);
%     y2(isub) = ycond(3) - ycond(2);
%     y2(isub) = ycond(3) - mean(ycond);
%     y2(isub) = (ycond(3) - ycond(2)) + (ycond(1) - ycond(2)) ;
%     y2(isub) = ycond(2) ;
%     y2(isub) = ycond(3) ;
%     y2(isub) = mean(ycond);
    predictor = [ (trt=='naloxone')*(-1) + (trt=='placebo')*0 + (trt=='morphine')*1 ];
    formula = ['y ~ 1 + opioid'];
    stat = fitglm(predictor,y,formula,'VarNames',{'opioid','y'});
    coef = stat.Coefficients;
    y2(isub) = coef.Estimate(2);

end

% correlation
clc;
[rho,p] = corr(y1,y2,'type','Spearman','rows','pairwise');
corrplot(y1,y2,'ctype','Spearman');
fprintf('Spearman correlation two-tailed test, rho = %.3g, p = %.2e \n',rho,p);
