%%  simulateBattery2
clear all;
clc;
close all;
          

%% set parameter space
mu = [-10:2:+10];
sigma = [0.001, 0.01 , 0.05, 0.1 , 0.5 , 1 , 5 , 10 , 50 , 100 , 1000];
alpha = [0.001, 0.01 , 0.05, 0.1 , 0.5 , 1 , 5 , 10 , 50 , 100 , 1000];
beta = [0.001, 0.01 , 0.05, 0.1 , 0.5 , 1 , 5 , 10 , 50 , 100 , 1000];

n_s = 1;
n_rep = 1;

for i_s = 1:n_s
      param(i_s).muV1 = 1; 
      param(i_s).sigmaV1 = 1;
      param(i_s).alpha = 1;
      param(i_s).beta = 2;
     
     param(i_s).muV2 = 1; param(i_s).sigmaV2 = 1;
     param(i_s).bR = 0;
     param(i_s).bE = 0;
     param(i_s).bm = 0;
     param(i_s).kRD = 30;
     param(i_s).kED = 30;
     param(i_s).eta = 10^4;

end

%% simulation loops

% initialization
medianR = nan(n_s,n_rep);
stdR = nan(n_s,n_rep);
medianE = nan(n_s,n_rep);
stdE = nan(n_s,n_rep);
predictedChoiceR = nan(n_s,n_rep);
predictedChoiceE = nan(n_s,n_rep);
acceptRE = nan(n_s,n_rep);
patientR = nan(n_s,n_rep);
patientE = nan(n_s,n_rep);


for i_s = 1:n_s
             
             
    % simulation function
    for j= 1:n_rep
        
        % sample hidden values
         param(i_s).V1 = safepos(param(i_s).sigmaV1*randn(1,24)+param(i_s).muV1);
         param(i_s).V2 = safepos(param(i_s).sigmaV2*randn(1,24)+param(i_s).muV2);
        
        % forward inference
        [y, u] = simBattery(param(i_s));

        % extract behavior
        ratingR = y(1,1:24);
        medianR(i_s,j) = median(ratingR);
        meanR(i_s,j) = mean(ratingR);
        stdR(i_s,j) = std(ratingR);

        
        ratingE = y(1,24+1:24+24);
        medianE(i_s,j) = median(ratingE);
        meanE(i_s,j) = mean(ratingE);
        stdE(i_s,j) = std(ratingE);

        choiceR = y(2,48+1:48+48);
        choiceE = y(2,96+1:96+48);
        dvR = param(i_s).V1(u(5,48+1:48+48)) - param(i_s).V1(u(4,48+1:48+48))   ;
        dvE = param(i_s).V2(u(5,96+1:96+48)) - param(i_s).V2(u(4,96+1:96+48))   ;
        predictedChoiceR(i_s,j) = mean( [mean(choiceR(dvR>0)) , 1- mean(choiceR(dvR<0))] );
        predictedChoiceE(i_s,j) = mean( [mean(choiceE(dvE<0)) , 1- mean(choiceE(dvE>0))] );
        
        
        weightRE = y(3,144+1:144+48);
        acceptRE(i_s,j) = mean(weightRE);
        
        discountR = y(4,192+1:192+48);
        patientR(i_s,j) = mean(discountR);

        
        discountE = y(4,240+1:240+48);
        patientE(i_s,j) = mean(discountE);
        
        % inverse inference
        [estimate(i_s,j)] = reverseBattery(y, u);
        



    end
end
             


%% display
% figure 1
     ind = [1:11];
     xlegend = ['mu_R'];
     p = [param(ind).muV1];
    
     f = figure; display_prediction_hiddenValuation;
     f = figure; display_inversion_hiddenValuation;

% figure 2
     ind = [12:22];
     xlegend = ['sigma_R'];
     p = [param(ind).sigmaV1];
    
     f = figure; display_prediction_hiddenValuation;
     f = figure; display_inversion_hiddenValuation;

    for i=1:numel(f.Children)
        f.Children(i).XScale = 'log';
    end
    
% figure 3
     ind = [23:33];
     xlegend = ['alpha'];
     p = [param(ind).alpha];
    
     f = figure; display_prediction_hiddenValuation;
     f = figure; display_inversion_hiddenValuation;

    for i=1:numel(f.Children)
        f.Children(i).XScale = 'log';
    end
    
% figure 4
     ind = [34:44];
     xlegend = ['beta'];
     p = [param(ind).beta];
    
     f = figure; display_prediction_hiddenValuation;
     f = figure; display_inversion_hiddenValuation;

    for i=1:numel(f.Children)
        f.Children(i).XScale = 'log';
    end



        
