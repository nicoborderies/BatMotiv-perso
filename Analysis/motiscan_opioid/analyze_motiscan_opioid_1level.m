%% analyze_motiscan_opioid_1level

% reset
clear all; close all; clc;

%% preparation
% folders definitions
datadir = 'B:\nicolas.borderies\projets\batmotiv\données\OPIOID';
codedir = 'B:\nicolas.borderies\projets\batmotiv\code.public\Bat-analyse';
subcodedir = 'B:\nicolas.borderies\projets\batmotiv\code.perso\Analysis\motiscan_opioid';
resultdir = 'B:\nicolas.borderies\projets\batmotiv\resultats\OPIOID';
addpath(genpath(codedir));
addpath(genpath(subcodedir));

% file definitions
resultname = 'motiscan-opioid-1level';

% load experimental design data
cd(datadir);
design = load('subject_table_opioid.mat');

% subject definition
sub4analysis = design.tab.n_inclusion(find(design.tab.selection==1))';
nsub = numel(sub4analysis);

%% option definition


%% analysis

% data preparation
groupdata = cell(1,nsub);
groupresult = cell(1,nsub);

for isub = sub4analysis % subject loop
    
    % display
    i = find(sub4analysis==isub);
    fprintf('analysis: subject %d / %d \n',i,nsub);
    
    % folder definition
    subdir = [datadir filesep 'sub' num2str(isub) ];
    cd(subdir);
     
    % load & format data
    [subdata] = load_gripAccu(subdir);
    
    % process data & get result
    [subresult,subdata] = process.gripAccu(subdata);
    
    % store
    groupdata{i} = subdata;
    groupresult{i} = subresult;
    
end

%% saving

cd(resultdir);
date = clock;
datename = [ num2str(date(3)) '_' num2str(date(2)) '_' num2str(date(1)) ];
save([resultname '-' datename '.mat'],'groupdata','groupresult','sub4analysis','design');




