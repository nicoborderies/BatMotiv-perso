%% analyze_batmotiv_opioid_1level

    clc;
    clear all;
    close all;


%% directory
    [root,vbadir,analysisdir,datadir,resultdir] = setPath;


%% specifications
    %subjects definition & selection
    setGroups;
    groups =  [ OPIOID  ];
    groupNames = {'OPIOID'};

%% options

    % individual options
    option = set_batmotiv_option;
    option.analysis.battery.set_analysis = @set_analysis_vRelease; % @set_analysis_vRelease;
    option.design.taskList =   {'gripAccu','grip','learning','rating','weight'};
%     option.design.taskList =   {'rating','weight'};

    option.analysis.learning.processor = @process_learning_opioid;
    option.analysis.gripAccu.selectedModel = 1;
    option.analysis.gripAccu.processor = @process_gripAccu;
    option.analysis.gripAccu.predictRest = 1;

    option.analysis.battery.predictRT=1;

    option.analysis.display=0;
    option.analysis.save=0;
    option.analysis.parallel=1;
    option.analysis.sequential=1;
    
%     % opioid
    if  ismember('OPIOID',groupNames)
        cd([datadir '\OPIOID']);
       set_design_opioid;
    end
    
    % subject loop options
    parallel = 0;
    compilation = 1;
    
    analysisName = 'batmotiv_analysis_OPIOID';



%% Execute individual analysis

% open parallel workers
    if parallel
        np = 7;
        parpool(np);
    end

    i=0; % initialize 
    % group loop
    for iGroup = 1:numel(groupNames)

        % define group data
            eval(['group = ' groupNames{iGroup} ';']);
            groupdir = [ datadir filesep groupNames{iGroup}];
            subject2analyze = group;
        % initialize
            cd([groupdir]);
            if exist( [groupNames{iGroup}],'file')==0;
                mkdir( [groupNames{iGroup}]) ;
            end % initialize

        % subject loop
        nSub = numel(subject2analyze);
        group_result  = cell(1,nSub);
        group_data  = cell(1,nSub);

        if parallel
            parfor iSub = 1:nSub
                   % display
                        fprintf('subject %d on %d \n',(i + iSub),numel(groups));

                   % analysis
                        subdir = [ groupdir filesep 'sub' num2str(subject2analyze(iSub))];
%                         try
                        [data,result] = analyze_batmotiv( subdir , option );
%                         end

                   % save & organize group data
                       group_result{iSub} = result;
                       group_data{iSub} = data;

                   % reset 
%                         clear  data result; % cannot clear because of transparency violation
                        close all; clc;

            end
        else
            for iSub = 1:nSub
                   % display
                        fprintf('subject %d on %d \n',(i + iSub),numel(groups));

                   % analysis
                        subdir = [ groupdir filesep 'sub' num2str(subject2analyze(iSub))];
                        [data,result] = analyze_batmotiv( subdir , option );

                   % save & organize group data
                       group_result{iSub} = result;
                       group_data{iSub} = data;

                   % reset 
                        clear  data result; close all; clc;
            end
        end

        % save & index independantly of the parallel processing
        for iSub = 1:nSub
           groupResult.(groupNames{iGroup}).subject{iSub} = group_result{iSub};
           groupData.(groupNames{iGroup}).subject{iSub} = group_data{iSub} ; 
           groupResult.population.subject{i + iSub} = group_result{iSub};
           groupData.population.subject{i + iSub} = group_data{iSub} ;
        end
        i = i + numel(subject2analyze);


    end
  
% close parallel workers
    if parallel
        delete(gcp);
    end
    
    
%% Compilation

if compilation
    % sub -> group
    for iG = 1:numel(groupNames)
        [ groupResult.(groupNames{iG}).group  , groupData.(groupNames{iG}).group  ] ...
            = compile_batmotiv_sub2group( groupResult.(groupNames{iG}).subject  , groupData.(groupNames{iG}).subject  , option );
    end

    % group -> population
    [ groupResult, groupData ] = compile_batmotiv_group2population( groupNames , groupResult, groupData , option  );
end

%% Save

% set path
    cd([resultdir]);
    
% save 
    date = clock;
    strDate = [ num2str(date(3)) '_' num2str(date(2)) '_' num2str(date(1)) ];
    
    
    save([ analysisName '_all_' strDate ],'groupResult','groupData','groupNames','groups','option');  
    
    if compilation
        dataStructure =  groupData.population  ;
        inferentialTable = groupResult.population.battery.inferential;
        misc = groupResult.population.battery.misc;
        save([ analysisName '_' strDate ],'dataStructure','inferentialTable','misc');  
    end
