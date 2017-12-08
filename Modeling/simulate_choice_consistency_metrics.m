%% simulate_choice_consistency_metrics

%% settings
% intialization
clear;
clc;
rng('shuffle');

% parameters
n_x = 24;
n_rep_item = 4;
nrep = n_x*n_rep_item;

n_iteration = 10;
beta_list = 1*10.^(-1:0.5:1);
n_rep_parameter = numel(beta_list);

% data preparation
report = struct;
report.random_accuracy = nan(1,n_rep_parameter*n_iteration);
report.consistent_accuracy = nan(1,n_rep_parameter*n_iteration);
report.random_entropy = nan(1,n_rep_parameter*n_iteration);
report.consistent_entropy = nan(1,n_rep_parameter*n_iteration);
report.beta = repelem(beta_list,n_iteration);

%% exploration of parameter space
icount = 0;
totalcount = n_rep_parameter*n_iteration;
tic;

for i_rep_param = 1:n_rep_parameter
    for iteration = 1:n_iteration 
        % counter
        icount = icount+1;
        t=toc; h = round(t/3600); min = round((t-h*3600)/60);
        clc;
        fprintf('# iteration = %d / %d \n',icount,totalcount);
        fprintf('# elapsed time = %d h %d min \n',h,min);

        % parameters
        beta = beta_list(i_rep_param);
        
        %%  simulations
        % true item ranks
        x1 = [];x2=[];
        for iRep = 1:n_rep_item
            x1 = [x1;randperm(n_x)'];
            x2 = [x2;randperm(n_x)'];
        end
        
        % - H0: random choices
        random_choices = double(randi(n_x,nrep,1) > n_x/2);

        % - H1: generative choices
        softmax = @(x) 1./(1+exp(-beta.*x));
        proba = softmax(x2-x1);
        consistent_choices = double(proba >= rand(nrep,1));

        %% consistency scores
        % H0: random choices
        % - choice-dependant ranking 
        item_value = find_rank(random_choices,x1,x2,n_x,nrep);
        % - entropy
        random_entropy = compute_entropy(item_value);
        % - classification accuracy
        random_accuracy = compute_accuracy(item_value,x1,x2,random_choices);


        % H1: generative choices
        % - choice-dependant ranking 
        item_value = find_rank(consistent_choices,x1,x2,n_x,nrep);
        % - entropy
        consistent_entropy = compute_entropy(item_value);
        % - classification accuracy
        consistent_accuracy = compute_accuracy(item_value,x1,x2,consistent_choices);
        
        %% store
        report.random_accuracy(iteration +(i_rep_param-1)*n_iteration) = random_accuracy;
        report.consistent_accuracy(iteration +(i_rep_param-1)*n_iteration) = consistent_accuracy;
        report.random_entropy(iteration +(i_rep_param-1)*n_iteration) = random_entropy;
        report.consistent_entropy(iteration +(i_rep_param-1)*n_iteration) = consistent_entropy;
    end
end

%% plots

% 1. cross-validation accuracy
f=figure; hold on; l =[];
l(1) = plot(splitapply(@nanmean,report.beta,findgroups(report.beta)),...
            splitapply(@nanmean,report.consistent_accuracy,findgroups(report.beta)),'b');
scatter(report.beta,report.consistent_accuracy,'filled','b');
l(2) = plot(splitapply(@nanmean,report.beta,findgroups(report.beta)),...
            splitapply(@nanmean,report.random_accuracy,findgroups(report.beta)),'Color',[1 1 1]*0.5);
scatter(report.beta,report.random_accuracy,'filled','MarkerFaceColor',[1 1 1]*0.5);
set_all_properties('LineWidth',2,'FontSize',12);
legend([l(1) l(2)],{'consistent choices','random choices'});
ylim([0 1]);
xlabel('beta (au.)');
ylabel('cross-classification accuracy (%)');

% 2. choice entropy
f=figure; hold on; l =[];
l(1) = plot(splitapply(@nanmean,report.beta,findgroups(report.beta)),...
            splitapply(@nanmean,report.consistent_entropy,findgroups(report.beta)),'b');
scatter(report.beta,report.consistent_entropy,'filled','b');
l(2) = plot(splitapply(@nanmean,report.beta,findgroups(report.beta)),...
            splitapply(@nanmean,report.random_entropy,findgroups(report.beta)),'Color',[1 1 1]*0.5);
scatter(report.beta,report.random_entropy,'filled','MarkerFaceColor',[1 1 1]*0.5);
set_all_properties('LineWidth',2,'FontSize',12);
legend([l(1) l(2)],{'consistent choices','random choices'});
ylim([0 n_x]);
xlabel('beta (au.)');
ylabel('choice entropy (bits)');




%% subfunctions

function item_value = find_rank(choice,x1,x2,n_x,nrep)

    chosen_item = (choice==0).*x1 + (choice==1).*x2;
    item_value = (histcounts(chosen_item,n_x)./histcounts([x1;x2],n_x))';

end

function entropy = compute_entropy(pd)
    pd(pd==0)=eps;
    pd(pd==1)=1-eps;
    entropy = - sum([pd.*log2(pd) + (1-pd).*log2(1-pd)]);
%     entropy = - mean([pd.*log2(pd) + (1-pd).*log2(1-pd)]);

end

function accuracy = compute_accuracy(item_value,x1,x2,choice)

    % circular method
%     predictor = item_value(x2)-item_value(x1);
%     subset = (predictor~=0);
%     accuracy = mean( (2*choice(subset)-1)==sign(predictor(subset)));
    
    % cross-validation method
    prediction = nan(numel(choice),1);
    for ichoice = 1:numel(choice)
        
        
        f1 = mean( [choice(setdiff(find(x1==x1(ichoice)),ichoice))==0 ;
                    choice(setdiff(find(x2==x1(ichoice)),ichoice))==1] );
        f2 = mean( [choice(setdiff(find(x1==x2(ichoice)),ichoice))==0 ;
                    choice(setdiff(find(x2==x2(ichoice)),ichoice))==1] );
        if sign(f2-f1)~=0
            prediction(ichoice) = double((f2-f1)>0);
        end
    end
    accuracy = mean(prediction(~isnan(prediction))==choice(~isnan(prediction)));
%     fprintf('accuracy = %f\n',accuracy);

end