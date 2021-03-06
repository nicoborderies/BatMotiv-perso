%% do_batmotiv_analysis_CONTROLS_sourceInformation

    clc;
    clear all;
    close all;


%% directory
    [root,vbadir,analysisdir,datadir,resultdir] = setPath;


%% specifications
    %subjects definition & selection
    setGroups;
    groups =  [ CONTROL  ];
    groupNames = {'CONTROL'};

%% options

    % individual options
    option = set_batmotiv_option;
    option.analysis.battery.set_analysis = @set_analysis_vRelease; % @set_analysis_vRelease;
    option.design.taskList =   {'rating','choice','weight','discount'};
    
%     option.analysis.battery.predict_rating=1;
%     option.analysis.battery.predict_choice=1;
%     option.analysis.battery.predict_weight=1;
%     option.analysis.battery.predict_discount=1;
%     option.analysis.battery.predict_rt=1;
    
    sourceDesign = struct2table(factorial_struct('rating',{0,1},'choice',{0,1},...
                                    'weight',{0,1},'discount',{0,1},'rt',{0,1}));
    sourceDesign = sourceDesign(2:end,:);
    source_rating = sourceDesign.rating;
    source_choice = sourceDesign.choice;
    source_rating = sourceDesign.weight;
    source_discount = sourceDesign.discount;
    source_rt = sourceDesign.rt;

    nsource = height(sourceDesign);
    
    option.analysis.display=0;
    option.analysis.save=0;
    option.analysis.parallel=1;
    option.analysis.sequential=0;
    
    
    % subject loop options
    parallel = 1;
    compilation = 0;
    
    analysisName = 'batmotiv_analysis_CONTROL_preferenceTasks_sourceInfo';



%% Execute individual analysis

% open parallel workers
    if parallel
        p = gcp('nocreate'); 
        if isempty(p)
            np = 7;
            parpool(np);
        end
        
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

       % citalopram
            if  ismember('CITALOPRAM',groupNames{iGroup})
               set_design_citalopram;
            else
                option.design.withinSubFactor = {'sessionNumber'};
                option.design.load_withinSubFactor = [];
            end

        % subject loop
        nSub = numel(subject2analyze);
%         group_result  = cell(nsource,nSub);
%         group_data  = cell(nsource,nSub);

        if parallel
            for is = 1:nsource
                    
               % option
                option.analysis.battery.predict_rating = sourceDesign.rating(is);
                option.analysis.battery.predict_choice = sourceDesign.choice(is);
                option.analysis.battery.predict_weight = sourceDesign.weight(is);
                option.analysis.battery.predict_discount = sourceDesign.discount(is);
                option.analysis.battery.predict_rt = sourceDesign.rt(is);
                    
                parfor iSub = 1:nSub
                   % display
                        fprintf('subject %d on %d \n',(i + iSub),numel(groups));

                   % analysis
                        subdir = [ groupdir filesep 'sub' num2str(subject2analyze(iSub))];
                        [data,result] = analyze_batmotiv( subdir , option );

                   % save & organize group data
                       group_result{is,iSub} = result;
                       group_data{is,iSub} = data;

                   % reset 
%                         clear  data result; % cannot clear because of transparency violation
                        close all; clc;
                end
            end
            
            
            
            
        else
            for iSub = 1:nSub
%             for iSub = jsub:nSub

               for is = 1:nsource
                    
                   % option
                    option.analysis.battery.predict_rating = sourceDesign.rating(is);
                    option.analysis.battery.predict_choice = sourceDesign.choice(is);
                    option.analysis.battery.predict_weight = sourceDesign.weight(is);
                    option.analysis.battery.predict_discount = sourceDesign.discount(is);
                    option.analysis.battery.predict_rt = sourceDesign.rt(is);
                    
                   % display
                        fprintf('subject %d on %d \n',(i + iSub),numel(groups));

                   % analysis
                        subdir = [ groupdir filesep 'sub' num2str(subject2analyze(iSub))];
                        [data,result] = analyze_batmotiv( subdir , option );

                   % save & organize group data
                       group_result{is,iSub} = result;
                       group_data{is,iSub} = data;

                   % reset 
%                         clear  data result; % cannot clear because of transparency violation
                        close all; clc;
                end

            end
        end

        % save & index independantly of the parallel processing
%         for iSub = 1:nSub
%            groupResult.(groupNames{iGroup}).subject{:,iSub} = group_result{:,iSub};
%            groupData.(groupNames{iGroup}).subject{:,iSub} = group_data{:,iSub} ; 
%            groupResult.population.subject{:,i + iSub} = group_result{:,iSub};
%            groupData.population.subject{:,i + iSub} = group_data{:,iSub} ;
%         end
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
    
    
    save([ analysisName '_all_' strDate ],'group_result','group_data','groupNames','groups','option');  
    save([ analysisName '_all_' strDate ],'group_result','group_data','groupNames','groups','option','-v7.3');  

    
    if compilation
        dataStructure =  groupData.population  ;
        inferentialTable = groupResult.population.battery.inferential;
        misc = groupResult.population.battery.misc;
        save([ analysisName '_' strDate ],'dataStructure','inferentialTable','misc');  
    end
