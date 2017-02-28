%% do_batmotiv_group_analysis

[root,vbadir,analysisdir,datadir,resultdir] = setPath;

% Load
    analysisName = 'batmotiv_controlPopulation_compare_14_1_2016.mat';
    load([ analysisName ]);
    cd(datadir);
    demographicTable = readtable('batmotivSubjectTable.xlsx');

    

%% Process Group


demographicTable.group = nominal(demographicTable.group);
demographicTable.sex = nominal(demographicTable.sex);
demographicTable.Properties.VariableNames{2} = 'subject';
demo = demographicTable(:,[1 2 5 6 7]);
inf = innerjoin(inferentialTable,demo);
% 

[designTable, statistics, correlations ] = process_batmotiv_group( inf, dataStructure, misc );


%% Visualise Group Result

displayCorrelations(correlations,designTable)

fig = struct;
for iTask = 1:numel(taskList) % iteration across tasks
    eval([ 'fig.' taskList{iTask} ,...
        ' = display_' opt.taskList{iTask} '_2(dataStructure.(' opt.taskList{iTask} ').table,[],inferentialTable,opt);' ]);
end



%% Save

% set path
    cd([resultdir]);

    save([ analysisName '_' strDate ],'demographicTable','inf','designTable','statistics','correlations','dimensionality','-append');  
