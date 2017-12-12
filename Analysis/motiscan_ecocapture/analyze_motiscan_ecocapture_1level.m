%% analyze_motiscan_ecocapture_1level

% set-up
clear all; 
close all; 
clc;
currentdir = pwd;

%% user-dependant parameters
% subid = input('subject number?\n');
subid = 6;


%% preparation
% file definitions
scriptname = 'analyze_motiscan_ecocapture_1level';
reportname = 'motiscan-ecocapture-1level';
reportname = [ reportname '-' datestr(clock,'dd-mmm-yyyy_HH.MM') '.mat' ];

% folders definitions
codedir = fileparts(which(scriptname));
datadir = 'B:\nicolas.borderies\owncloud\motiscan\motiscan-ecocapture\motiscan-ecocapture-data';
reportdir = 'B:\nicolas.borderies\projets\batmotiv\reports\ecocapture';
reportdir = [reportdir filesep scriptname];
if exist(reportdir)~=7; mkdir(reportdir);end

% subject-level anaylsis definitions
taskList = {'ratingR','ratingE','choiceR','choiceE','weighRE',...
            'choiceMoneyDelay','choiceMoneyProba',...
            'rating_grip','ConfidencePrecision','gripRP'};
        
% graphical parameters
set(groot,'defaultLineLineWidth',1.5);
set(groot,'defaultAxesFontSize',12);

        
% file detection
subdir = [ datadir filesep 'sub' num2str(subid) ];
cd(subdir);
f = dir;
if numel(f)<=2; warning('file detection: not enough data files!');end

%% analysis

subdata = struct;

%% 1. Preference tests
% -------------------
% 1.1. Ratings

% file detector
find_testname  = @(dim) ['*rating' dim '*sub' num2str(subid) '*'];

%% 1.1.1. Rewards
% -------------------------
% parameters
dim = 'R';
dimname = 'reward';

% file detection/extraction
filepattern = find_testname(dim);
testfile = dir(filepattern);
d = load(testfile.name);

% formating
% - enumeration
ntrain = size(d.training_data,1);
ntest = size(d.data,1);
ntrial = ntrain + ntest;
% - variables creation 
phase = [ repmat(nominal('train'),ntrain,1) ; repmat(nominal('test'),ntest,1) ];
dimension = repmat(nominal(dimname),ntrain+ntest,1);
sessionNumber = [d.training_data(:,1);d.data(:,1)];
trialNumber = [d.training_data(:,2);d.data(:,2)];
itemIndex = [d.training_data(:,2);d.data(:,3)];
itemName = [nominal(d.traininglist(d.training_data(:,2))');nominal(d.([dimname 'list'])(d.data(:,3))')];
rating_discrete = round(( [d.training_data(:,3);d.data(:,4)] )./10)./10;
rating_continuous = ( [d.training_data(:,3);d.data(:,4)] )./100;
response_time = [d.training_data(:,4);d.data(:,5)];
confirmation_time = [d.training_data(:,5);d.data(:,6)];
% - table creation
subdata.rating.datatab = table(phase,dimension,sessionNumber,trialNumber,itemIndex,itemName,...
                               rating_discrete,rating_continuous,response_time,confirmation_time);
                           
% display
f=figure;
nquantile = 11;
histogram(subdata.rating.datatab.rating_discrete(subdata.rating.datatab.dimension==dimname),nquantile,...
          'Normalization','Probability');
xlim([0 1]);
xlabel('rating (%)');
ylabel('frequency (%)');
title(dimname);
     

%% 1.1.2. Efforts
% -------------------------
% parameters
dim = 'E';
dimname = 'effort';

% file detection/extraction
filepattern = find_testname(dim);
testfile = dir(filepattern);
d = load(testfile.name);

% formating
% - enumeration
ntrain = size(d.training_data,1);
ntest = size(d.data,1);
ntrial = ntrain + ntest;
% - variables creation 
phase = [ repmat(nominal('train'),ntrain,1) ; repmat(nominal('test'),ntest,1) ];
dimension = repmat(nominal(dimname),ntrain+ntest,1);
sessionNumber = [d.training_data(:,1);d.data(:,1)];
trialNumber = [d.training_data(:,2);d.data(:,2)];
itemIndex = [d.training_data(:,2);d.data(:,3)];
itemName = [nominal(d.traininglist(d.training_data(:,2))');nominal(d.([dimname 'list'])(d.data(:,3))')];
rating_discrete = round(( [d.training_data(:,3);d.data(:,4)] )./10)./10;
rating_continuous = ( [d.training_data(:,3);d.data(:,4)] )./100;
response_time = [d.training_data(:,4);d.data(:,5)];
confirmation_time = [d.training_data(:,5);d.data(:,6)];
% - table creation
subdata.rating.datatab =  [subdata.rating.datatab;
                           table(phase,dimension,sessionNumber,trialNumber,itemIndex,itemName,...
                               rating_discrete,rating_continuous,response_time,confirmation_time)];

% display
f=figure;
nquantile = 11;
histogram(subdata.rating.datatab.rating_discrete(subdata.rating.datatab.dimension==dimname),nquantile,...
          'Normalization','Probability');
xlim([0 1]);
xlabel('rating (%)');
ylabel('frequency (%)');
title(dimname);


% statistics
subset = (subdata.rating.datatab.phase=='test');
substat.rating.descriptive = varfun(@nanmean,subdata.rating.datatab(subset,:),'InputVariables',{'rating_continuous'},'GroupingVariables',{'dimension'});
substat.rating.descriptive = join(substat.rating.descriptive,...
                             varfun(@nanvar,subdata.rating.datatab(subset,:),'InputVariables',{'rating_continuous'},'GroupingVariables',{'dimension'}));
substat.rating.descriptive = join(substat.rating.descriptive,...
                             varfun(@entropy,subdata.rating.datatab(subset,:),'InputVariables',{'rating_discrete'},'GroupingVariables',{'dimension'}));
print_my_table(substat.rating.descriptive);


%% 1.2. Subjective choices

% file detector
find_testname  = @(dim) ['*choice' dim '*sub' num2str(subid) '*'];

%% 1.1.1. Rewards
% ----------------------
% parameters
dim = 'R';
dimname = 'reward';

% file detection/extraction
filepattern = find_testname(dim);
testfile = dir(filepattern);
d = load(testfile.name);


% formating
% - enumeration
ntrain = size(d.training_data,1);
ntest = size(d.data,1);
ntrial = ntrain + ntest;
% - variables creation 
phase = [ repmat(nominal('train'),ntrain,1) ; repmat(nominal('test'),ntest,1) ];
dimension = repmat(nominal(dimname),ntrain+ntest,1);
sessionNumber = [d.training_data(:,1);d.data(:,1)];
trialNumber = [d.training_data(:,2);d.data(:,2)];
leftItemIndex = [[1 3 5 1 2 3]';d.data(:,4)];
rightItemIndex = [[2 4 6 4 5 6]';d.data(:,6)];
leftRating_continuous = [nan(ntrain,1);d.data(:,5)./100];
rightRating_continuous = [nan(ntrain,1);d.data(:,7)./100];
sideChoice = [d.training_data(:,3);d.data(:,8)];
% sideChoicePosition = [d.training_data(:,4);d.data(:,9)]; % NOTE: to include for future datasets
response_time = [d.training_data(:,4);d.data(:,9)];
% - table creation
subdata.choice.datatab =  [table(phase,dimension,sessionNumber,trialNumber,leftItemIndex,rightItemIndex,leftRating_continuous,rightRating_continuous,...
                               sideChoice,response_time)];

% display
% - prepare graphic objects
f=figure;
nquantile = 6;
% - data conditionalization
subset = (subdata.choice.datatab.phase=='test' & subdata.choice.datatab.dimension==dimname);
dv = subdata.choice.datatab.rightRating_continuous(subset)-subdata.choice.datatab.leftRating_continuous(subset);
choice = double(subdata.choice.datatab.sideChoice(subset)==1);
rt = (subdata.choice.datatab.response_time(subset));
dv_bin = quantileranks(dv,nquantile);
% - choice curve
subplot(1,2,1);
x_mean = splitapply(@nanmean,dv,dv_bin);
y_mean = splitapply(@nanmean,choice,dv_bin);
y_sem = splitapply(@sem,choice,dv_bin);
[~,~,h] = errorscat( x_mean  , y_mean, y_sem,'k');
xlabel('\delta rating (right-left) (%)');
ylabel('choice = right (%)');
title(dimname);
% - response time curve
subplot(1,2,2);
x_mean = splitapply(@nanmean,dv,dv_bin);
y_mean = splitapply(@nanmean,rt,dv_bin);
y_sem = splitapply(@sem,rt,dv_bin);
[~,~,h] = errorscat( x_mean  , y_mean, y_sem,'k');
xlabel('\delta rating (right-left) (%)');
ylabel('response time (s.)');


% statistics
% - accuracy score
accuracy = mean([ mean(choice(dv>0)) , mean(1-choice(dv<0)) ]);
fprintf('Accuracy (%s) = %f \n',dim,accuracy);

%% 1.1.2. Efforts
% ----------------------
% parameters
dim = 'E';
dimname = 'effort';

% file detection/extraction
filepattern = find_testname(dim);
testfile = dir(filepattern);
d = load(testfile.name);

% formating
% - enumeration
ntrain = size(d.training_data,1);
ntest = size(d.data,1);
ntrial = ntrain + ntest;
% - variables creation 
phase = [ repmat(nominal('train'),ntrain,1) ; repmat(nominal('test'),ntest,1) ];
dimension = repmat(nominal(dimname),ntrain+ntest,1);
sessionNumber = [d.training_data(:,1);d.data(:,1)];
trialNumber = [d.training_data(:,2);d.data(:,2)];
leftItemIndex = [[1 3 5 1 2 3]';d.data(:,4)];
rightItemIndex = [[2 4 6 4 5 6]';d.data(:,6)];
leftRating_continuous = [nan(ntrain,1);d.data(:,5)./100];
rightRating_continuous = [nan(ntrain,1);d.data(:,7)./100];
sideChoice = [d.training_data(:,3);d.data(:,8)];
% sideChoicePosition = [d.training_data(:,4);d.data(:,9)]; % NOTE: to include for future datasets
response_time = [d.training_data(:,4);d.data(:,9)];
% - table creation
subdata.choice.datatab =  [subdata.choice.datatab;
                           table(phase,dimension,sessionNumber,trialNumber,leftItemIndex,rightItemIndex,leftRating_continuous,rightRating_continuous,...
                               sideChoice,response_time)];

% display
% - prepare graphic objects
f=figure;
nquantile = 6;
% - data conditionalization
subset = (subdata.choice.datatab.phase=='test' & subdata.choice.datatab.dimension==dimname);
dv = subdata.choice.datatab.rightRating_continuous(subset)-subdata.choice.datatab.leftRating_continuous(subset);
choice = double(subdata.choice.datatab.sideChoice(subset)==1);
rt = (subdata.choice.datatab.response_time(subset));
dv_bin = quantileranks(dv,nquantile);
% - choice curve
subplot(1,2,1);
x_mean = splitapply(@nanmean,dv,dv_bin);
y_mean = splitapply(@nanmean,choice,dv_bin);
y_sem = splitapply(@sem,choice,dv_bin);
[~,~,h] = errorscat( x_mean  , y_mean, y_sem,'k');
xlabel('\delta rating (right-left) (%)');
ylabel('choice = right (%)');
title(dimname);
% - response time curve
subplot(1,2,2);
x_mean = splitapply(@nanmean,dv,dv_bin);
y_mean = splitapply(@nanmean,rt,dv_bin);
y_sem = splitapply(@sem,rt,dv_bin);
[~,~,h] = errorscat( x_mean  , y_mean, y_sem,'k');
xlabel('\delta rating (right-left) (%)');
ylabel('response time (s.)');


% statistics
% - accuracy score
accuracy = mean([ mean(1-choice(dv>0)) , mean(choice(dv<0)) ]);
fprintf('Accuracy (%s) = %f \n',dim,accuracy);



%% 1.1.3. Rewards/Efforts
% ----------------------------------

% file detector
find_testname  = @() ['*weightRE*sub' num2str(subid) '*'];

% parameters
benefitname = 'reward';
costname = 'effort';

% file detection/extraction
filepattern = find_testname();
testfile = dir(filepattern);
d = load(testfile.name);

% formating
% - enumeration
ntrain = size(d.training_data,1);
ntest = size(d.data,1);
ntrial = ntrain + ntest;
% - variables creation 
phase = [ repmat(nominal('train'),ntrain,1) ; repmat(nominal('test'),ntest,1) ];
benefit = repmat(nominal(benefitname),ntrain+ntest,1);
cost = repmat(nominal(costname),ntrain+ntest,1);
sessionNumber = [d.training_data(:,1);d.data(:,1)];
trialNumber = [d.training_data(:,2);d.data(:,2)];
benefitItemIndex = [[1 2 3 4 5 6]';d.data(:,3)];
costItemIndex = [[1 2 3 4 5 6]';d.data(:,5)];
benefitRating_continuous = [nan(ntrain,1);d.data(:,4)./100];
costRating_continuous = [nan(ntrain,1);d.data(:,6)./100];
sideChoice = [d.training_data(:,3);d.data(:,7)];
% sideChoicePosition = [d.training_data(:,4);d.data(:,9)]; % NOTE: to include for future datasets
acceptChoice = [d.training_data(:,4);d.data(:,8)];
response_time = [d.training_data(:,5);d.data(:,9)];
% - table creation
subdata.choice2.datatab =  [table(phase,benefit,cost,sessionNumber,trialNumber,benefitItemIndex,costItemIndex,benefitRating_continuous,costRating_continuous,...
                               sideChoice,acceptChoice,response_time)];

% display
% - prepare graphic objects
f=figure;
nquantile = 6;
% - data conditionalization
subset = (subdata.choice2.datatab.phase=='test');
benefit = subdata.choice2.datatab.benefitRating_continuous(subset);
cost = subdata.choice2.datatab.costRating_continuous(subset);
choice = double(subdata.choice2.datatab.acceptChoice(subset)==1);
rt = (subdata.choice2.datatab.response_time(subset));
% - choice curves
% -- benefit elasticity
subplot(1,2,1);
x_bin = quantileranks(benefit,nquantile);
x_mean = splitapply(@nanmean,benefit,x_bin);
y_mean = splitapply(@nanmean,choice,x_bin);
y_sem = splitapply(@sem,choice,x_bin);
[~,~,h] = errorscat( x_mean  , y_mean, y_sem,'k');
xlabel('benefit rating  (%)');
ylabel('choice = right (%)');
title(benefitname);
% -- cost elasticity
subplot(1,2,2);
x_bin = quantileranks(cost,nquantile);
x_mean = splitapply(@nanmean,cost,x_bin);
y_mean = splitapply(@nanmean,choice,x_bin);
y_sem = splitapply(@sem,choice,x_bin);
[~,~,h] = errorscat( x_mean  , y_mean, y_sem,'k');
xlabel('cost rating  (%)');
ylabel('choice = right (%)');
title(costname);
                           
%% 1.3. Objective choices

% file detector
find_testname  = @(dim) ['*choiceMoney' dim '*sub' num2str(subid) '*'];

%% 1.3.1. Money/Delay
% ----------------------

% parameters
dim = 'Delay';
dimname = 'delay';

% file detection/extraction
filepattern = find_testname(dim);
testfile = dir(filepattern);
d = load(testfile.name);

% formating
% - enumeration
ntrain = size(d.training_data,1);
ntest = size(d.data,1);
ntrial = ntrain + ntest;
% - variables creation 
phase = [ repmat(nominal('train'),ntrain,1) ; repmat(nominal('test'),ntest,1) ];
dimension = repmat(nominal(dimname),ntrain+ntest,1);
sessionNumber = [d.training_data(:,1);d.data(:,1)];
trialNumber = [d.training_data(:,2);d.data(:,2)];
easyAmount = [d.training_data(:,3);d.data(:,3)];
hardAmount = [d.training_data(:,4);d.data(:,4)];
easyCost_ordinal = [d.training_data(:,5);d.data(:,5)];
hardCost_ordinal = [d.training_data(:,7);d.data(:,7)];
easyCost = [d.training_data(:,6);d.data(:,6)];
hardCost = [d.training_data(:,8);d.data(:,8)];
dimensionality = [d.training_data(:,9);d.data(:,9)];
obvious = double(easyAmount>hardAmount);
sideHard = [d.training_data(:,10);d.data(:,10)];
hardChoice = double([d.training_data(:,12);d.data(:,11)]==1); % NOTE: inversion to modify??
sideChoice = [d.training_data(:,11);d.data(:,12)];
% sideChoicePosition = [d.training_data(:,4);d.data(:,9)]; % NOTE: to include for future datasets
response_time =  [d.training_data(:,13);d.data(:,13)];
% - table creation
subdata.choice3.datatab =  [table(phase,dimension,sessionNumber,trialNumber,...
                                easyAmount,hardAmount,easyCost_ordinal,hardCost_ordinal,easyCost,hardCost,dimensionality,obvious,sideHard,...
                               sideChoice,hardChoice,response_time)];

% display
% - prepare graphic objects
f=figure;
nquantile = 6;
% - data conditionalization
subset = (subdata.choice3.datatab.phase=='test' & subdata.choice3.datatab.dimension==dimname & subdata.choice3.datatab.obvious==0);
distractor = subdata.choice3.datatab.easyAmount(subset);
cost = subdata.choice3.datatab.hardCost(subset);
easy_amount = subdata.choice3.datatab.easyAmount(subset);
hard_amount = subdata.choice3.datatab.hardAmount(subset);
choice = double(subdata.choice3.datatab.hardChoice(subset)==1);
chosen_amount = hard_amount;
chosen_amount(choice==0) = easy_amount(choice==0);
dimensionality = subdata.choice3.datatab.dimensionality(subset);
% - choice curves
% -- temporal indifference
for idim = 1:2
    x_bin = findgroups(cost(dimensionality==idim));
    x_mean = splitapply(@nanmean,cost(dimensionality==idim),x_bin);
    y_mean = splitapply(@nanmean,chosen_amount(dimensionality==idim),x_bin);
    y_sem = splitapply(@sem,chosen_amount(dimensionality==idim),x_bin);
    [~,~,h(idim)] = errorscat( x_mean  , y_mean, y_sem,'k');
    if idim==2
        h(idim).LineStyle='--'; 
    end
end
% ylim([0 50]);
legend([h(1) h(2)],'1-D','2-D')
xlabel('cost (delay)');
ylabel('equivalent money (€)');
title([ dimname '-indifference curve' ]);

                           
%% 1.3.2. Money/Probability

% parameters
dim = 'Proba';
dimname = 'probability';

% file detection/extraction
filepattern = find_testname(dim);
testfile = dir(filepattern);
d = load(testfile.name);

% formating
% - enumeration
ntrain = size(d.training_data,1);
ntest = size(d.data,1);
ntrial = ntrain + ntest;
% - variables creation 
phase = [ repmat(nominal('train'),ntrain,1) ; repmat(nominal('test'),ntest,1) ];
dimension = repmat(nominal(dimname),ntrain+ntest,1);
sessionNumber = [d.training_data(:,1);d.data(:,1)];
trialNumber = [d.training_data(:,2);d.data(:,2)];
easyAmount = [d.training_data(:,3);d.data(:,3)];
hardAmount = [d.training_data(:,4);d.data(:,4)];
easyCost_ordinal = [d.training_data(:,5);d.data(:,5)];
hardCost_ordinal = [d.training_data(:,7);d.data(:,7)];
easyCost = [d.training_data(:,6);d.data(:,6)];
hardCost = [d.training_data(:,8);d.data(:,8)];
dimensionality = [d.training_data(:,9);d.data(:,9)];
obvious = double(easyAmount>hardAmount);
sideHard = [d.training_data(:,10);d.data(:,10)];
hardChoice = double([d.training_data(:,12);d.data(:,11)]==1); % NOTE: inversion to modify??
sideChoice = [d.training_data(:,11);d.data(:,12)];
% sideChoicePosition = [d.training_data(:,4);d.data(:,9)]; % NOTE: to include for future datasets
response_time =  [d.training_data(:,13);d.data(:,13)];
% - table creation
subdata.choice3.datatab =  [subdata.choice3.datatab;
                            table(phase,dimension,sessionNumber,trialNumber,...
                            easyAmount,hardAmount,easyCost_ordinal,hardCost_ordinal,easyCost,hardCost,dimensionality,obvious,sideHard,...
                            sideChoice,hardChoice,response_time)];

% display
% - prepare graphic objects
f=figure;
nquantile = 6;
% - data conditionalization
subset = (subdata.choice3.datatab.phase=='test' & subdata.choice3.datatab.dimension==dimname & subdata.choice3.datatab.obvious==0);
distractor = subdata.choice3.datatab.easyAmount(subset);
cost = subdata.choice3.datatab.hardCost(subset);
easy_amount = subdata.choice3.datatab.easyAmount(subset);
hard_amount = subdata.choice3.datatab.hardAmount(subset);
choice = double(subdata.choice3.datatab.hardChoice(subset)==1);
chosen_amount = hard_amount;
chosen_amount(choice==0) = easy_amount(choice==0);
dimensionality = subdata.choice3.datatab.dimensionality(subset);
% - choice curves
% -- cost indifference
for idim = 1:2
    x_bin = findgroups(cost(dimensionality==idim));
    x_mean = splitapply(@nanmean,cost(dimensionality==idim),x_bin);
    y_mean = splitapply(@nanmean,chosen_amount(dimensionality==idim),x_bin);
    y_sem = splitapply(@sem,chosen_amount(dimensionality==idim),x_bin);
    [~,~,h(idim)] = errorscat( x_mean  , y_mean, y_sem,'k');
    if idim==2
        h(idim).LineStyle='--'; 
    end
end
% ylim([0 50]);
legend([h(1) h(2)],'1-D','2-D')
xlabel('cost (probability)');
ylabel('equivalent money (€)');
title([ dimname '-indifference curve' ]);


%% 2. Force tests 
%------------------------
%% 2.1. Force perception

% file detector
find_testname  = @() ['*rating_grip*sub' num2str(subid) '*'];

% file detection/extraction
filepattern = find_testname();
testfile = dir(filepattern);
d = load(testfile.name);

% formating
% - enumeration
ntest = size(d.data,1);
ntrial = ntest;
% - variables creation 
phase = [ repmat(nominal('test'),ntest,1) ];
sessionNumber = [d.data(:,1)];
trialNumber = [d.data(:,2)];
forceLevel = [d.data(:,3)];
force = [d.data(:,4)];
force_normalized = [d.data(:,5)./100];
rating = [d.data(:,6)./100];
response_time = [d.data(:,7)];
confirmation_time = [d.data(:,8)];
% - table creation
subdata.forcePerception.datatab =  [table(phase,sessionNumber,trialNumber,forceLevel,force,force_normalized,rating,...
                                    response_time,confirmation_time)];

% display
% - prepare graphic objects
f=figure;
% - data conditionalization
instruction = subdata.forcePerception.datatab.forceLevel;
force = subdata.forcePerception.datatab.force_normalized;
perception = subdata.forcePerception.datatab.rating;
% - raw data curves
% -- force control
subplot(1,2,1);hold on;
[ h,hp,s ] = corrplot(instruction,force);
% axis([0 1 0 1]);
xlabel('target force (%fmax)');
ylabel('exerted force (%fmax)');
title([ 'force control' ]);
% -- force perception
subplot(1,2,2);hold on;
[ h,hp,s ] = corrplot(force,perception);
% axis([0 1 0 1]);
xlabel('exerted force (%fmax)');
ylabel('perceived penibility (%)');
title([ 'effort perception' ]);
                                
%% 2.2. Force control

% file detector
find_testname  = @() ['*ConfidencePrecision*sub' num2str(subid) '*'];

% file detection/extraction
filepattern = find_testname();
testfile = dir(filepattern);
d = load(testfile.name);

% formating
% - enumeration
% - enumeration
ntrain = size(d.training_data,1);
ntest = size(d.data,1);
ntrial = ntrain + ntest;
% - variables creation 
datatab=struct;
% -- temporal structure
datatab.phase = [ repmat(nominal('train'),ntrain,1) ; repmat(nominal('test'),ntest,1) ];
datatab.sessionNumber = [ ones(ntrain,1) ; ones(ntest,1)].*d.nsession ;
datatab.trialNumber = [d.training_data(:,1);d.data(:,1)];
datatab.blockNumber = [nan(ntrain,1);d.data(:,2)];
% -- conditions 
datatab.feedbackBias = [zeros(ntrain,1);d.data(:,4)];
datatab.feedbackValidity = [ones(ntrain,1);d.data(:,5)];
datatab.incentive = [nan(ntrain,1);d.data(:,3)];
% -- performances
datatab.force = [d.training_data(:,9);d.data(:,7)];
datatab.error = [d.training_data(:,10);d.data(:,8)];
datatab.absolute_error = [d.training_data(:,11);d.data(:,9)];
datatab.correctForce = [d.training_data(:,12);d.data(:,10)];
datatab.correctTime = [d.training_data(:,13);d.data(:,11)];
datatab.correct =  [ double(d.training_data(:,12) & d.training_data(:,13)) ; d.data(:,12)];
datatab.feedback = [d.training_data(:,4);d.data(:,6)];
datatab.gain =  [d.training_data(:,3);d.data(:,13)];
% -- subjective estimates
datatab.confidence = [nan(ntrain,1);d.data(:,14)]./100;
datatab.confidence_confirmation_time = [nan(ntrain,1);d.data(:,16)];
datatab.mood = [nan(ntrain,1);d.data(:,17)]./100;
datatab.mood_confirmation_time = [nan(ntrain,1);d.data(:,19)];
% - table creation
subdata.forceControl.datatab =  struct2table(datatab);


% display
% - prepare graphic objects
f=figure;
% - data conditionalization
subset = (subdata.forceControl.datatab.phase=='test');
feedbackBias = subdata.forceControl.datatab.feedbackBias(subset);
incentive = subdata.forceControl.datatab.incentive(subset);
force = subdata.forceControl.datatab.error(subset);
confidence = subdata.forceControl.datatab.confidence(subset);
% - raw data curves
% -- force control
xvar = [feedbackBias,incentive,feedbackBias,incentive];
yvar = [force,force,confidence,confidence];
xvarlegend = {'feedback condition','incitation (€)','feedback condition','incitation (€)'};
yvarlegend = {'force error (%fmax)','force error (%fmax)','confidence (%)','confidence (%)'};
for iplot=1:4
    subplot(2,2,iplot);hold on;
    x_bin = findgroups(xvar(:,iplot));
    x_mean = splitapply(@nanmean,xvar(:,iplot),x_bin);
    y_mean = splitapply(@nanmean,yvar(:,iplot),x_bin);
    y_sem = splitapply(@sem,yvar(:,iplot),x_bin);
    [~,~,h] = errorscat( x_mean  , y_mean, y_sem,'k');
    xlabel(xvarlegend{iplot});
    ylabel(yvarlegend{iplot});
    set(gca,'XTick',x_mean);
    if mod(iplot,2)==0
       set(gca,'XScale','log') 
    end
end


%% 2.3. Force power

% file detector
find_testname  = @() ['*gripRP*sub' num2str(subid) '*'];

% file detection/extraction
filepattern = find_testname();
testfile = dir(filepattern);
d = load(testfile.name);

% formating
% - enumeration
% - enumeration
ntrain = size(d.training_data,1);
ntest = size(d.data,1);
ntrial = ntrain + ntest;
% - variables creation 
datatab=struct;
% -- temporal structure
datatab.phase = [ repmat(nominal('train'),ntrain,1) ; repmat(nominal('test'),ntest,1) ];
datatab.sessionNumber = [ ones(ntrain,1) ; ones(ntest,1)].*d.nsession ;
datatab.trialNumber = [d.training_data(:,1);d.data(:,1)];
% -- conditions 
datatab.incentive = [d.training_data(:,2);d.data(:,2)];
% -- performances
datatab.force = [d.training_data(:,3);d.data(:,3)];
datatab.forceCumulative = [d.training_data(:,4);d.data(:,4)];
datatab.force_normalized = [d.training_data(:,5);d.data(:,5)]./100;
datatab.gain =  [d.training_data(:,6);d.data(:,6)];
% - table creation
subdata.forcePower.datatab =  struct2table(datatab);

% display
% - prepare graphic objects
f=figure;
nquantile=10;
% - data conditionalization
subset = (subdata.forcePower.datatab.phase=='test');
trialNumber = subdata.forcePower.datatab.trialNumber(subset);
incentive = subdata.forcePower.datatab.incentive(subset);
force = subdata.forcePower.datatab.force_normalized(subset);
% - force curves
% -- incentive modulation
subplot(1,2,1);hold on;
x_bin = findgroups(incentive);
x_mean = splitapply(@nanmean,incentive,x_bin);
y_mean = splitapply(@nanmean,force,x_bin);
y_sem = splitapply(@sem,force,x_bin);
[~,~,h] = errorscat( x_mean  , y_mean, y_sem,'k');
xlabel('incentive level');
ylabel('exerted force (%fmax)');
set(gca,'XTick',x_mean);
% -- fatigue modulation
subplot(1,2,2);hold on;
x = quantileranks(trialNumber,nquantile);
x_bin = findgroups(x);
x_mean = round(splitapply(@nanmean,trialNumber,x_bin));
y_mean = splitapply(@nanmean,force,x_bin);
y_sem = splitapply(@sem,force,x_bin);
[~,~,h] = errorscat( x_mean  , y_mean, y_sem,'k');
xlabel('trial index');
ylabel('exerted force (%fmax)');
set(gca,'XTick',x_mean);

%% save 

cd(reportdir);
save(reportname,'subdata','-mat');
% makewordreport( scriptname,reportdir);
% writetable(subdata.rating.datatable,'ratingtable.xlsx');

%% reset
% path reset
cd(currentdir);