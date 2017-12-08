%% norris_test
%
% Author:  Nicolas Borderies
% email address: nico.borderies@gmail.com 
% May 2017; Last revision: 


clear all;
clc;

% identification
nsub = input([ 'subject inclusion number?\n' ]);
nsess = input([ 'session number?\n' ]);

% load questions and answers
testdir = 'B:\nicolas.borderies\projets\batmotiv\code.perso\Analysis\motiscan_opioid';
datadir = 'B:\nicolas.borderies\projets\batmotiv\données\OPIOID';
subdir = [ datadir '\sub' num2str(nsub) ];
cd(testdir);
items = load('t_Norris.mat');

% define variables
ntrial = size(items.t_Norris,1);
rating = nan(ntrial,1);
norm = 9.55;

% parameters
conf_instruct = 'En ce moment, je me sens:\n';

% run trials
clc; disp(['Appuyer sur une touche pour débuter le test ']); pause;

for t=1:ntrial
    
    % EVA report
    bounds = {items.t_Norris{t,:}};
    fprintf(['\n']);
    disp([' ']);
    rating(t) = input([  bounds{1} ' / '  bounds{2} ' ?\n' ]);
    clc;
     
end

% compute scores
score = rating/norm;
mean_score = nanmean(score);


% save 
cd(subdir);
boundnames = items.t_Norris;
filename = [ 'norris_score' '_sub' num2str(nsub) '_sess' num2str(nsess)];
save(filename,'rating','boundnames','score','mean_score');
cd(testdir);


