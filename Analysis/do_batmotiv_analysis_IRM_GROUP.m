%% do_batmotiv_analysis_IRM_GROUP

    clc;
    clear all;
    close all;


%% directory
    [root,vbadir,analysisdir,datadir,resultdir] = setPath;


%% specifications
    %subjects definition & selection
    setGroups;
    groups =  [ IRM  ];
    groupNames = {'IRM'};

%% options

    % individual options
    option = set_batmotiv_option;
    option.analysis.battery.set_analysis = @set_analysis_vRelease; % @set_analysis_vRelease;
    option.design.taskList =   {'rating','choice','weight'};
    option.design.group = 'IRM';

    option.analysis.display=0;
    option.analysis.save=0;
    option.analysis.parallel=0;
    option.analysis.sequential=0;

%     % citalopram
%     if  ismember('CITALOPRAM',groupNames)
%        set_design_citalopram;
%     end
    
    % subject loop options
    parallel = 1;
    
    analysisName = 'batmotiv_analysis_IRM_preferenceTasks';



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
            
            subFileList = dir; 
            fileNameList = {subFileList(3:end).name};

       % citalopram
            if  ismember('CITALOPRAM',groupNames{iGroup})
               set_design_citalopram;
            else
                option.design.withinSubFactor = {'sessionNumber'};
                option.design.load_withinSubFactor = [];
            end

        % subject loop
        nSub = numel(subject2analyze);
        group_result  = cell(1,nSub);
        group_data  = cell(1,nSub);

        if parallel
            parfor iSub = 1:nSub
                   % display
                        fprintf('subject %d on %d \n',(i + iSub),numel(groups));

                   % analysis
%                         subdir = [ groupdir filesep 'sub' num2str(subject2analyze(iSub))];
                        subdir = [ groupdir filesep fileNameList{iSub} filesep 'behavior' ];
                        [data,result] = analyze_batmotiv_irm( subdir , option );

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
%                         subdir = [ groupdir filesep 'sub' num2str(subject2analyze(iSub))];
                        subdir = [ groupdir filesep fileNameList{iSub} filesep 'behavior' ];
                        [data,result] = analyze_batmotiv_irm( subdir , option );

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

% sub -> group
for iG = 1:numel(groupNames)
    [ groupResult.(groupNames{iG}).group  , groupData.(groupNames{iG}).group  ] ...
        = compile_batmotiv_sub2group( groupResult.(groupNames{iG}).subject  , groupData.(groupNames{iG}).subject  , option );
end

% group -> population
[ groupResult, groupData ] = compile_batmotiv_group2population( groupNames , groupResult, groupData , option  );


%% Save

% set path
    cd([resultdir]);
    
% save 
    date = clock;
    strDate = [ num2str(date(3)) '_' num2str(date(2)) '_' num2str(date(1)) ];
    
    
    save([ analysisName '_all_' strDate ],'groupResult','groupData','groupNames','groups','option');  
    
    dataStructure =  groupData.population  ;
    inferentialTable = groupResult.population.battery.inferential;
    misc = groupResult.population.battery.misc;
    save([ analysisName '_' strDate ],'dataStructure','inferentialTable','misc');  
