function  [data,result,option,fig]=analyze_batmotiv_irm( subdir , option )
% analyze_batmotiv 
%
% this functions execute an individual analysis of the complete batmotiv
% battery as follows: 
%
%   -INPUTS:
%       - subdir: a string indicating the directory path where the
%                   subject data is located
%       - option: a structure specifying the method of analysis
%               with the following fields
%           - design
%           - analysis
%
%
%   -OUTPUTS:
%       - data: a structure storing raw & preprocessed data, 
%               with the following fields
%           - battery
%           - 'taskName' (one field for every element of option.design.taskList ) 
%           
%       - result: a structure storing descriptive & inferential statistics
%                 estimated based on the data, details of the model used 
%                 with the following fields
%           - battery
%           - 'taskName' (one field for every element of option.design.taskList ) 
%
%
%       - fig: a structure storing figures visualizing the analyzing scheme
%                 with the following fields
%           - 'taskName' (one field for every element of option.design.taskList ) 
%
%
%
% Nicolas Borderies
% 11/2015
%
%
%% Default option
cd(subdir);
if nargin<2
      [option] = set_batmotiv_option;
else
      [option] = set_batmotiv_option(option);
end

%% Specifications
cd(subdir);

    %  design specification
        taskList = option.design.taskList;
        taskSerial = {'grip','gripIAPS','mental','learning'};
        
        
    % analysis specification
        set_analysis = option.analysis.battery.set_analysis;
        [ ~ , option] = set_analysis(option);

%% Initialization 
cd(subdir);

    % directories
        ind = strfind(subdir,'sub');
        subid =  str2num(subdir(ind+4:ind+7));

    % mat file
        if exist( ['batmotiv_sub' num2str(subid) '.mat'])==2;
            if option.analysis.reload==0
                % reset current analysis
                delete( ['batmotiv_sub' num2str(subid) '.mat']) ; 
                result=struct; result.battery =struct;  data=struct;
            else
                % import previous analysis
                load( ['batmotiv_sub' num2str(subid) '.mat']) ; 
            end
        else
            result=struct;result.battery=struct;  data=struct;
        end 

        
%% Data Analysis

% set path
cd(subdir);


    % iteration across tasks
    for iTask = 1:numel(taskList) 
        
                fprintf('analyze task:  %s  \n',taskList{iTask});

                % Load
                if ~isfield(data,taskList{iTask})
                    loader = str2func(['load_' taskList{iTask} '_irm' ]);
                    [data.(taskList{iTask})] = loader(subdir,option);
                end

                % Process sequentially task analysis
                if option.analysis.sequential || (option.analysis.parallel && ismember(taskList{iTask},taskSerial))
                    if ~isempty(find(data.(taskList{iTask}).misc.nSession~=0))
                        processor = str2func(['process_' taskList{iTask} ]);
                        [result.(taskList{iTask}),data.(taskList{iTask}),option] = processor(data.(taskList{iTask}),option);
                    end
                    
                    % Text
                    if isfield(result,taskList{iTask})
                        if isfield(result.(taskList{iTask}),'inferential')
                           fprintf('\n');
                           disp(table2struct(result.(taskList{iTask}).inferential))
                        end
                    end
                
                end
                
    end
    
    

    % whole-battery analysis
        % compile
        [result,data,option] = compile_battery(subid,result,data,option);

        % process whole-battery
        if option.analysis.parallel
            fprintf('analyze battery  \n');

            [result.battery , data.battery , option] = process_battery( data.battery  , option );
            
            % update predictions
                for iTask = 1:numel(taskList) 
                    data.(taskList{iTask}).table = data.battery.table(data.battery.table.task==taskList{iTask},:);
                end
        end
                        
         % Text
               fprintf('model list: \n');
               disp(option.analysis.battery.modelList)

            if isfield(result.battery,'inferential')
               fprintf('results: \n');
               disp(table2struct(result.battery.inferential))
            end
        
            
            
            
    % Display
        for iTask = 1:numel(taskList) 
            if option.analysis.display
                displayer = str2func(['display_' taskList{iTask} ]);
                dd = data.(taskList{iTask});
                try
                   rr = result.(taskList{iTask});
                catch
                   rr = []; 
                end
                [fig,result.(taskList{iTask}),data.(taskList{iTask})] = displayer(rr,dd,option);
            end
        end


%% Saving Analysis
    % set path
clc;
cd(subdir);

    % matfiles
    if option.analysis.save
        save(['batmotiv_sub' num2str(subid)],'data','result','option');
    end


end

