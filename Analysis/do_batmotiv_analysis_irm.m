%% do_batmotiv_analysis_irm

    clc;
    clear all;
    close all;


%% directory
    [root,vbadir,analysisdir,datadir,resultdir] = setPath;


%% specifications
    %subjects definition & selection
    setGroups;
    subdir = [datadir '\IRM\sub_8001_250416\behavior'] ;

%% options
    option = set_batmotiv_option;
    option.analysis.battery.set_analysis = @set_analysis_vRelease; % @set_analysis_vRelease;
    option.design.taskList =   {'rating','choice','weight'};
    option.design.group = 'IRM';

    option.analysis.display=0;
    option.analysis.save=0;
    option.analysis.parallel=0;
    option.analysis.sequential=0;
    
    % citalopram
    if ~isempty(strfind(subdir,'CITALOPRAM'))
       set_design_citalopram;
    end


%% Execute individual analysis

    [data,result,option] = analyze_batmotiv_irm( subdir , option );

     