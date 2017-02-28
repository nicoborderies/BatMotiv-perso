%%  simulateBattery3
clear all;
clc;
close all;
          

%% set parameter space
n_s = 12;
n_rep = 10;
nt = 48;

mu = ones(1,n_s*n_rep);
sigma = ones(1,n_s*n_rep);
beta = ones(1,n_s*n_rep);
eta = 1e4.*ones(1,n_s*n_rep);

mu(1:30) = repmat([ 0.5  1  2 ],1,n_rep);
sigma(31:60) = repmat([ 0.5  1  2 ],1,n_rep);
beta(61:90) = repmat([ 0.5  1  2 ],1,n_rep);
eta(91:120) = repmat([ 1e3 , 1e4 , 1e5 ],1,n_rep);


%% SECTION TITLE
% DESCRIPTIVE TEXT


for i_s = 1:n_s
        for j= 1:n_rep

         ind = j + (i_s-1)*n_rep ;  

         % generative parameters
         param(ind).muV1 = mu(i_s); 
         param(ind).sigmaV1 = sigma(i_s);
         param(ind).alpha = 1;
         param(ind).beta = beta(i_s);

         param(ind).muV2 =  0.3;
         param(ind).sigmaV2 = 1;

         param(ind).bR = 0;
         param(ind).bE = 0;
         param(ind).bm = 0;
         param(ind).kRD = 1;
         param(ind).kED = 30;
         param(ind).eta = eta(i_s);
         
         % predictions
        prediction(ind).medianR = NaN;
        prediction(ind).stdR = NaN;
        prediction(ind).medianE = NaN;
        prediction(ind).stdE = NaN;
        prediction(ind).predictedChoiceR = NaN;
        prediction(ind).predictedChoiceE = NaN;
        prediction(ind).acceptRE = NaN;
        prediction(ind).patientR = NaN;
        prediction(ind).patientE = NaN;
        prediction(ind).discountR = nan(1,nt);
        prediction(ind).delay = nan(1,nt);
         
        end
end

clear mu sigma alpha beta ind;

%% simulation loops

for i_s = 1:n_s
             
    % simulation function
    for j= 1:n_rep
        
         ind = j + (i_s-1)*n_rep ;  

        
        % sample hidden values with respect to a parametric law
%             % normal 
%              param(ind).V1 = (param(ind).sigmaV1.*randn(1,24) + param(ind).muV1);
%              param(ind).V2 = (param(ind).sigmaV2.*randn(1,24) + param(ind).muV2);           
            
%             % log-normal 
%              param(ind).V1 = exp(param(ind).sigmaV1.*randn(1,24) + param(ind).muV1);
%              param(ind).V2 = exp(param(ind).sigmaV2.*randn(1,24) + param(ind).muV2);
         
            % safepos
%              param(ind).V1 = safepos(param(ind).sigmaV1.*randn(1,24) + param(ind).muV1);
%              param(ind).V2 = safepos(param(ind).sigmaV2.*randn(1,24) + param(ind).muV2);
             
%             % pos-normal 
             param(ind).V1 = sqrt((param(ind).sigmaV1.*randn(1,24) + param(ind).muV1).^2);
             param(ind).V2 = sqrt((param(ind).sigmaV2.*randn(1,24) + param(ind).muV2).^2);

%             % gamma
%              pd1 = makedist('Gamma','a',param(ind).muV1,'b',param(ind).sigmaV1);
%              pd2 = makedist('Gamma','a',param(ind).muV2,'b',param(ind).sigmaV2);
%              param(ind).V1 = random(pd1,[1 24]);
%              param(ind).V2 = random(pd2,[1 24]);
           
            
        
        % forward inference
        [y, u] = simBattery(param(ind));

        % extract behavior
            % rating
            ratingR = y(1,1:24);
            ratingE = y(1,24+1:24+24);

            prediction(ind).medianR = median(ratingR);
            prediction(ind).meanR = mean(ratingR);
            prediction(ind).stdR = std(ratingR);
            prediction(ind).medianE = median(ratingE);
            prediction(ind).meanE = mean(ratingE);
            prediction(ind).stdE = std(ratingE);

            % choice 1D
            choiceR = y(2,48+1:48+48);
            choiceE = y(2,96+1:96+48);
            dvR = param(ind).V1(u(5,48+1:48+48)) - param(ind).V1(u(4,48+1:48+48))   ;
            dvE = param(ind).V2(u(5,96+1:96+48)) - param(ind).V2(u(4,96+1:96+48))   ;
            
            prediction(ind).predictedChoiceR = mean( [mean(choiceR(dvR>0)) , 1- mean(choiceR(dvR<0))] );
            prediction(ind).predictedChoiceE = mean( [mean(choiceE(dvE<0)) , 1- mean(choiceE(dvE>0))] );
        
            % choice 2D
            weightRE = y(3,144+1:144+48);
            
            prediction(ind).acceptRE = mean(weightRE);

            % choice IT 
            discountR = y(4,192+1:192+48);
            discountE = y(4,240+1:240+48);

            prediction(ind).discountR(:) = y(4,192+1:192+48);
            prediction(ind).delay(:) = u(17,192+1:192+48);
            prediction(ind).patientR= mean(discountR);
            prediction(ind).patientE = mean(discountE);
        
        % inverse inference
        try
           [estimate(ind)] = reverseBattery(y, u);
        end
        
        % clear
        clear y u ratingR ratingE choiceR choiceE dvR dvE weightRE discountR discountE ;


    end
end
             
%% saving
save('simulations_HiddenValue2','param','prediction');
try
    save('simulations_HiddenValue2','estimate','-append') ;
end

%% display
% 
    f = figure; hold on;
    
    subplot(1,4,1);
    for i = 1:numel(param)
        x(i) = std(param(i).V1)./mean(param(i).V1);
        dev = abs(param(i).V1-median(param(i).V1));
        x(i) =  mean(dev)./mean(param(i).V1) ;
        x(i) =  madiff(param(i).V1,0) ;
        
        y(i) = prediction(i).predictedChoiceR ;
    end
    [r,p] = corr(x',y','type','Pearson');
    scatter(x,y,30,'MarkerFaceColor','k','MarkerEdgeColor','k');
    l = lsline;
    text(0.8*max(x),1.1*min(y),['r = ' num2str(r)],'Color','r','FontWeight','bold');
    text(0.8*max(x),1.05*min(y),['p = ' num2str(p)],'Color','r','FontWeight','bold');

    l.Color = 'r';
    xlabel('mean normalized distance ');
    ylabel('predictable choice (1D)');

    subplot(1,4,2);
    for i = 1:numel(param)
        x(i) = iqr(param(i).V1);
        y(i) = prediction(i).predictedChoiceR ;
    end
    [r,p] = corr(x',y','type','Pearson');
    scatter(x,y,30,'MarkerFaceColor','k','MarkerEdgeColor','k');
    l = lsline;
    text(0.8*max(x),1.1*min(y),['r = ' num2str(r)],'Color','r','FontWeight','bold');
    text(0.8*max(x),1.05*min(y),['p = ' num2str(p)],'Color','r','FontWeight','bold');
    l.Color = 'r';
    xlabel('range');
    ylabel('predictable choice (1D)');
    
    subplot(1,4,3);
    for i = 1:numel(param)
        x(i) = kurtosis(param(i).V1);
        y(i) = prediction(i).predictedChoiceR ;
    end
    [r,p] = corr(x',y','type','Pearson');
    scatter(x,y,30,'MarkerFaceColor','k','MarkerEdgeColor','k');
    l = lsline;
    text(0.8*max(x),1.1*min(y),['r = ' num2str(r)],'Color','r','FontWeight','bold');
    text(0.8*max(x),1.05*min(y),['p = ' num2str(p)],'Color','r','FontWeight','bold');
    l.Color = 'r';
    xlabel('kurtosis');
    ylabel('predictable choice (1D)');
    
    subplot(1,4,4);
    for i = 1:numel(param)
%         x(i) = entropy_shannon(param(i).V1,0.01,[0 10]);
        x(i) = param(i).beta ;

        y(i) = prediction(i).predictedChoiceR ;
    end
    [r,p] = corr(x',y','type','Pearson');
    scatter(x,y,30,'MarkerFaceColor','k','MarkerEdgeColor','k');
    l = lsline;
    text(0.8*max(x),1.1*min(y),['r = ' num2str(r)],'Color','r','FontWeight','bold');
    text(0.8*max(x),1.05*min(y),['p = ' num2str(p)],'Color','r','FontWeight','bold');
    l.Color = 'r';
    xlabel('entropy');
    ylabel('predictable choice (1D)');
    
    
    
    f = figure; hold on;
    
    subplot(1,3,1);
    for i = 1:numel(param)
        x(i) = median(param(i).V1);
        y(i) = prediction(i).acceptRE ;
    end
    [r,p] = corr(x',y','type','Pearson');
    scatter(x,y,30,'MarkerFaceColor','k','MarkerEdgeColor','k');
    l = lsline;
    text(0.8*max(x),1.1*min(y),['r = ' num2str(r)],'Color','r','FontWeight','bold');
    text(0.8*max(x),1.05*min(y),['p = ' num2str(p)],'Color','r','FontWeight','bold');

    l.Color = 'r';
    xlabel('median ');
    ylabel('accept choice (2D)');
    
    subplot(1,3,2);
    for i = 1:numel(param)
        x(i) = median(param(i).V1)^2/(madiff(param(i).V1)^2);
        y(i) = prediction(i).acceptRE ;
    end
    [r,p] = corr(x',y','type','Pearson');
    scatter(x,y,30,'MarkerFaceColor','k','MarkerEdgeColor','k');
    l = lsline;
    text(0.8*max(x),1.1*min(y),['r = ' num2str(r)],'Color','r','FontWeight','bold');
    text(0.8*max(x),1.05*min(y),['p = ' num2str(p)],'Color','r','FontWeight','bold');

    l.Color = 'r';
    xlabel('rate');
    ylabel('accept choice (2D)');
    
    
    subplot(1,3,3);
    for i = 1:numel(param)
        x(i) = median(param(i).V1) -  median(param(i).V2);
        y(i) = prediction(i).acceptRE ;
    end
    [r,p] = corr(x',y','type','Pearson');
    scatter(x,y,30,'MarkerFaceColor','k','MarkerEdgeColor','k');
    l = lsline;
    text(0.8*max(x),1.1*min(y),['r = ' num2str(r)],'Color','r','FontWeight','bold');
    text(0.8*max(x),1.05*min(y),['p = ' num2str(p)],'Color','r','FontWeight','bold');

    l.Color = 'r';
    xlabel('mode');
    ylabel('accept choice (2D)');
    
    
% % figure predictions
%     d  = mean(delay,2);
%     choice = mean(discountR,2);
%     
%     f = figure; hold on;
%     x = tools.tapply(delay,{delay},@nanmean,{'continuous'},8);
%     y = tools.tapply(discountR,{delay},@nanmean,{'continuous'},8);
%     z = tools.tapply(discountR,{delay},@sem,{'continuous'},8);
%     
%     h = errorbar(x,y,z,'LineStyle','-','Color','k','Marker','o') ;
%     xlabel('delays');
%     ylabel('patient choice (%)');
% 
% % % figure 1
% %      ind = [1:3];
% %      xlegend = ['mu_R'];
% %      p = [param(ind).muV1];
% %     
% %      f = figure; display_prediction_hiddenValuation;
% %      f = figure; display_inversion_hiddenValuation;
% 
% % figure 2
%      ind = [1:3];
%      xlegend = ['sigma_R'];
%      p = [param(ind).sigmaV1];
%     
%      f = figure; display_prediction_hiddenValuation;
% %      f = figure; display_inversion_hiddenValuation;
% 
%     for i=1:numel(f.Children)
%         f.Children(i).XScale = 'log';
%     end
% %     
% % % figure 3
% %      ind = [23:33];
% %      xlegend = ['alpha'];
% %      p = [param(ind).alpha];
% %     
% %      f = figure; display_prediction_hiddenValuation;
% %      f = figure; display_inversion_hiddenValuation;
% % 
% %     for i=1:numel(f.Children)
% %         f.Children(i).XScale = 'log';
% %     end
% %     
% % % figure 4
% %      ind = [34:44];
% %      xlegend = ['beta'];
% %      p = [param(ind).beta];
% %     
% %      f = figure; display_prediction_hiddenValuation;
% %      f = figure; display_inversion_hiddenValuation;
% % 
% %     for i=1:numel(f.Children)
% %         f.Children(i).XScale = 'log';
% %     end



        
