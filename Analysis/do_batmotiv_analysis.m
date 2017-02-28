%% demo_batmotiv_analysis
% this script demonstrate a subject-level batmotiv analysis scheme
%


clc;
clear all;
close all;


% pathways set-up
%%% (personalize wrt. your path names)
root = [ 'B\nicolas.borderies\Projets scientifiques\Projets BatMotiv'];
datadir  = [ root filesep 'Bat-données'];
subdir = [ datadir '\CONTROL\sub13'] ;

% options
option = set_batmotiv_option; % this function set the default options
option.analysis.battery.set_analysis = @set_analysis_vRelease; % @set_analysis_vRelease;
option.design.taskList =   {'rating','choice','weight','discount'}; % list of task to analyze
option.analysis.display=0; % flag to display task figures
option.analysis.save=0; % flag to save the results
option.analysis.parallel=1;  % flag to perform a parallel/aggregated analysis of multiple tasks 
                             %(until now, only preference tasks under the hidden valuation model)
option.analysis.sequential=0;  % flag to perform a sequential/segregated analysis of multiple tasks 

% run analysis
[data,result,option] = analyze_batmotiv( subdir , option );

     