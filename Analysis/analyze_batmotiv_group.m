%% analyze_batmotiv_group

clc;
clear all;
close all;


%% Directory Configuration

    [root,vbadir,analysisdir,datadir,resultdir] = setPath;


%% Specifications

    %subjects definition & selection
    setGroups;
            

%% Options

% opt.taskList = {'rating','choice','weight','discount','learning','grip','gripIAPS','mental'};
opt.taskList = {'rating','choice','weight','discount'};

opt.reload = 0; % flag to indicate wheter to reset the analysis (0) or reload (1) 
opt.parallel = 1; % parallel processing


%% Execute individual analysis

if opt.parallel
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

        if opt.parallel
            parfor iSub = 1:nSub
                   % display
                        fprintf('subject %d on %d \n',(i + iSub),numel(groups));

                   % analysis
                        subdir = [ groupdir filesep 'sub' num2str(subject2analyze(iSub))];
                        if opt.reload==0
                            % reset
                            [data,result]=analyze_batmotiv( subdir );
                        elseif opt.reload==1
                            % reload
                            file = load( [ subdir '\batmotiv_sub' num2str(subject2analyze(iSub)) '.mat']) ;
                            data = file.data;  result = file.result;
                        end

                   % save & organize group data
                       group_result{iSub} = result;
                       group_data{iSub} = data;
                       
                   % reset 
%                         clear  data result; 
                        close all; clc;

            end
        else
            for iSub = 1:nSub
                   % display
                        fprintf('subject %d on %d \n',(i + iSub),numel(groups));

                   % analysis
                        subdir = [ groupdir filesep 'sub' num2str(subject2analyze(iSub))];
                        if opt.reload==0
                            % reset
                            [data,result]=analyze_batmotiv( subdir );
                        elseif opt.reload==1
                            % reload
                            load( [ subdir '\batmotiv_sub' num2str(subject2analyze(iSub)) '.mat']) ;
                        end

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
        
    if opt.parallel
        delete(gcp);
    end
        
%% Compilation

% sub -> group
for iG = 1:numel(groupNames)
    [ groupResult.(groupNames{iG}).group  , groupData.(groupNames{iG}).group  ] ...
        = compile_batmotiv_sub2group( groupResult.(groupNames{iG}).subject  , groupData.(groupNames{iG}).subject  , opt );
end

% group -> population
[ groupResult, groupData ] = compile_batmotiv_group2population( groupNames , groupResult, groupData , opt  );


   


%% Save

% set path
    cd([resultdir]);
    
% save 
    analysisName = 'batmotiv_control_ftd_hiddenValue';
    date = clock;
    strDate = [ num2str(date(3)) '_' num2str(date(2)) '_' num2str(date(1)) ];
    
    
    save([ analysisName '_all_' strDate ],'groupResult','groupData','groupNames','groups','opt');  
    
    dataStructure =  groupData.population  ;
    inferentialTable = groupResult.population.battery.inferential;
    misc = groupResult.population.battery.misc;
    save([ analysisName '_' strDate ],'dataStructure','inferentialTable','misc');  





     