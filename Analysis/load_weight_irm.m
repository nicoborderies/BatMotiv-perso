
function [data] = load_weight_irm(subjectDir,option)
%% INITIALIZE 
    % variable definition
        data=struct;
        if nargin>1;
           if isequal( class(option.design.load_withinSubFactor) , 'function_handle' )
               load_withinSubFactor = option.design.load_withinSubFactor;
           else
               load_withinSubFactor = [];
           end
        else
            load_withinSubFactor = [];
        end

    
    % task selection
        subManipList = {'WeightRE_s', 'WeightPE_s', 'WeightRP_s'}; % to use as MANIP_NAME
        benefitItemList={'rewardlist', 'punishlist', 'rewardlist'};            % to get list of item
        costItemList={'effortlist', 'effortlist', 'punishlist'};            % to get list of item

        benefitFieldName={'writtenReward', 'writtenPunishment', 'writtenReward'}; % to use as field name in data
        costFieldName={'writtenEffort', 'writtenEffort', 'writtenPunishment'}; % to use as field name in data
        benefit = {'reward', 'punishment', 'reward'}; % to use as field name in data
        cost = {'effort', 'effort', 'punishment'}; % to use as field name in data
        
        submanipFieldName={'RewardEffort', 'PunishmentEffort', 'RewardPunishment'}; % to use as field name in data
        submanipAcronym = {'RE','PE','RP'};
        nSubManip=length(submanipFieldName);
        data.misc = table(zeros(nSubManip,1),'RowNames',submanipFieldName,'VariableNames',{'nSession'});

    % session selection
        sessionList = [2:4]; 
        
%% ITERATION ACROSS TASKS
for iSubmanip = 1:nSubManip
    %% FILE EXTRACTION
    % ===========================================================================
    MANIP_NAME = subManipList{iSubmanip}; % result files should include manip 'manipName'
    
    % test for file accessibility
     fileList=dir([subjectDir filesep '*' MANIP_NAME '*.mat']);
    if isempty(fileList)
        warning(['No file found for ' MANIP_NAME 'in directory: ' subjectDir])
        isFileFound =0;
    else
        isFileFound =1;
    end
    
        
    if isFileFound == 1 % check file existence
        [fileList, nSession] = tools.load_file(subjectDir, MANIP_NAME,sessionList);
        data.misc{submanipFieldName{iSubmanip},'nSession'} = nSession;
        data.misc{submanipFieldName{iSubmanip},'acronym'} = {submanipAcronym{iSubmanip}};
        
        % concatenate sessions
        allData = vertcat(fileList.data);
        
        %% MISC
        % ===========================================================================
        
        % --- correction for erroneous saving of rating*Dim*summary --- %
        get_correctRating_weight;
        data.misc{submanipFieldName{iSubmanip},'fileName'} = {fileList.fileName};
        data.misc{submanipFieldName{iSubmanip},'benefitList'} = {fileList.([benefitItemList{iSubmanip}])};
        data.misc{submanipFieldName{iSubmanip},'costList'} = {fileList.([costItemList{iSubmanip}])};

        % rating summaries
         switch submanipFieldName{iSubmanip}
            case {'RewardEffort'}
                  data.misc{submanipFieldName{iSubmanip},'benefitRating'} = {fileList.ratingRsummary};
                  data.misc{submanipFieldName{iSubmanip},'costRating'} = {fileList.ratingEsummary};
            case {'PunishmentEffort'}
                  data.misc{submanipFieldName{iSubmanip},'benefitRating'} = {fileList.ratingPsummary};
                  data.misc{submanipFieldName{iSubmanip},'costRating'} = {fileList.ratingEsummary};
            case {'RewardPunishment'}
                  data.misc{submanipFieldName{iSubmanip},'benefitRating'} = {fileList.ratingRsummary};
                  data.misc{submanipFieldName{iSubmanip},'costRating'} = {fileList.ratingPsummary};
         end
         

        %% BEHAVIOUR
        % ===========================================================================
        
        subdata.behavior.RT  = allData(:, 6)';                         % RT
        subdata.behavior.sideChoice  = allData(:, 7)';                 % choice side: [1] right [-1] left
        subdata.behavior.isGoChoice  = allData(:, 8)';                 % choice type: [1] go : accept the trade-off [0] no go : refuse the trade-off
        
        
        %% CONDITIONS
        % ===========================================================================
        subdata.condition.sessionNumber=[];
        for iSession = 1:nSession
            [nTrial,~]=size(fileList(iSession).data);
            subdata.condition.sessionNumber =   [subdata.condition.sessionNumber iSession * ones(1,nTrial)] ;
        end
        
        
        
        subdata.condition.trialNumber                = allData(:, 1)';
        subdata.condition.benefitItemNumber          = allData(:, 2)';                         % benefit item number : reward in RP & RE, (avoid) punishment in PE
        subdata.condition.costItemNumber             = allData(:, 4)';                         % cost item number : effort in RE & PE, punishment in RP
        subdata.condition.benefitItemRatingValue     = allData(:, 3)';                         % rating of the benefit item number : reward in RP & RE, (avoid) punishment in PE
        subdata.condition.costItemRatingValue        = allData(:, 5)';                         % rating of the cost item number : reward in RP & RE, (avoid) punishment in PE
        
        subdata.condition.sideGo = 1*(((subdata.behavior.sideChoice == 1) & (subdata.behavior.isGoChoice == 1)) | ((subdata.behavior.sideChoice == -1) & (subdata.behavior.isGoChoice == 0)));
        subdata.condition.sideGo( subdata.condition.sideGo~=1)= -1;
        
%         itemsName = {'benefitItemName','costItemName'};
%         benefitList = data.misc.benefitList{1,1}; 
%         costList = data.misc.costList{1,1}; 
%         [l,~] = size(benefitList);
%         if l~=24; benefitList = benefitList';end     
%         [l,~] = size(costList);
%         if l~=24; costList = costList';end  
%         subdata.condition.benefitItemName =  nominal(benefitList(subdata.condition.benefitItemNumber)');
%         subdata.condition.costItemName =  nominal(costList(subdata.condition.costItemNumber)');
% 
%         
%         %  variable transformations
%          b_list  = data.misc{submanipFieldName{iSubmanip},'benefitRating'}{:};
%          c_list = data.misc{submanipFieldName{iSubmanip},'costRating'}{:};
%          muBenefit = mean(b_list);
%          sigmaBenefit=std(b_list);
%          muCost=mean(c_list);
%          sigmaCost=std(c_list);
%          
%          subdata.condition.benefitStandardValue = (subdata.condition.benefitItemRatingValue)/sigmaBenefit;
%          subdata.condition.costStandardValue  = (subdata.condition.costItemRatingValue)/sigmaCost; 
%          subdata.condition.benefitItemZValue = (subdata.condition.benefitItemRatingValue-muBenefit)/sigmaBenefit;
%          subdata.condition.costItemZValue  = (subdata.condition.costItemRatingValue-muCost)/sigmaCost; 
%          subdata.condition.benefitItemRangeValue = (subdata.condition.benefitItemRatingValue-min(b_list))/range(b_list);
%          subdata.condition.costItemRangeValue  = (subdata.condition.costItemRatingValue-min(c_list))/range(c_list);
         
         
         % consensual ratings
%          try
%             cd ..
%             cd ..
%             load('consensualValue.mat');
%             cd(subjectDir);
%             [~,ind] = ismember(subdata.condition.benefitItemName',consensus.itemName);
%             subdata.condition.benefitConsensusValue = consensus.value(ind)';
%             [~,ind] = ismember(subdata.condition.costItemName',consensus.itemName);
%             subdata.condition.costConsensusValue = consensus.value(ind)';
%             subdata.condition.rewardPersonnalValue = zeros(1,numel(subdata.condition.benefitConsensusValue));
%             subdata.condition.punishmentPersonnalValue = zeros(1,numel(subdata.condition.benefitConsensusValue));
%             subdata.condition.effortPersonnalValue = zeros(1,numel(subdata.condition.benefitConsensusValue));
%             
%             % search sub index
%                 groupName = option.design.group;
%                 subind = strfind(subjectDir,'sub');
%                 subid = str2num((subjectDir(subind+3:end)));
%                 ind = find(scaling.group==groupName &...
%                            scaling.subject==subid &...
%                            scaling.dimension==costFieldName{iSubmanip})  ;
%   
%                 gain = scaling.gain(ind);
%                 offset = scaling.offset(ind);
%                 subdata.condition.costPersonnalValue = (subdata.condition.costItemRatingValue./100)./gain - offset; 
%                 subdata.condition.([ cost{iSubmanip} 'PersonnalValue' ]) = subdata.condition.costPersonnalValue;
%         end
         
         % dimensions
%          switch submanipFieldName{iSubmanip}
%             case {'RewardEffort'}
%                   subdata.condition.rewardRatingValue = subdata.condition.benefitItemRatingValue;
%                   subdata.condition.punishmentRatingValue = zeros(1,numel(subdata.condition.costItemZValue));
%                   subdata.condition.effortRatingValue = subdata.condition.costItemRatingValue;
%                   
%                   subdata.condition.rewardItemZValue = subdata.condition.benefitItemZValue;
%                   subdata.condition.punishmentItemZValue = zeros(1,numel(subdata.condition.costItemZValue));
%                   subdata.condition.effortItemZValue = subdata.condition.costItemZValue;
%                   
%                   subdata.condition.rewardStandardValue = subdata.condition.benefitStandardValue;
%                   subdata.condition.punishmentStandardValue = zeros(1,numel(subdata.condition.costItemZValue));
%                   subdata.condition.effortStandardValue = subdata.condition.costStandardValue;
%              case {'PunishmentEffort'}
%                   subdata.condition.rewardRatingValue = zeros(1,numel(subdata.condition.costItemRatingValue));
%                   subdata.condition.punishmentRatingValue = subdata.condition.benefitItemRatingValue;
%                   subdata.condition.effortRatingValue = subdata.condition.costItemRatingValue;
%                   
%                   subdata.condition.rewardItemZValue = zeros(1,numel(subdata.condition.costItemZValue));
%                   subdata.condition.punishmentItemZValue = subdata.condition.benefitItemZValue;
%                   subdata.condition.effortItemZValue = subdata.condition.costItemZValue;
%                   
%                   subdata.condition.rewardStandardValue = zeros(1,numel(subdata.condition.costItemZValue));
%                   subdata.condition.punishmentStandardValue = subdata.condition.benefitStandardValue;
%                   subdata.condition.effortStandardValue = subdata.condition.costStandardValue;
%             case {'RewardPunishment'}
%                   subdata.condition.rewardRatingValue = subdata.condition.benefitItemRatingValue;
%                   subdata.condition.punishmentRatingValue = subdata.condition.costItemRatingValue;
%                   subdata.condition.effortRatingValue = zeros(1,numel(subdata.condition.costItemZValue));
%                   
%                   subdata.condition.rewardItemZValue = subdata.condition.benefitItemZValue;
%                   subdata.condition.punishmentItemZValue = subdata.condition.costItemZValue;
%                   subdata.condition.effortItemZValue = zeros(1,numel(subdata.condition.costItemZValue));
%                   
%                   subdata.condition.rewardStandardValue = subdata.condition.benefitStandardValue;
%                   subdata.condition.punishmentStandardValue = subdata.condition.costStandardValue;
%                   subdata.condition.effortStandardValue = zeros(1,numel(subdata.condition.costItemZValue));
%          end
%         
          % within-subject factor
        if ~isempty(load_withinSubFactor)
            [subdata.condition] =  load_withinSubFactor(subdata.condition,subjectDir,option);
        end
        
        % design quality control
        [rho, p ] = corr(subdata.condition.benefitItemRatingValue',subdata.condition.costItemRatingValue');
        if p<=0.05; h=1; else h=0;end
        data.misc{submanipFieldName{iSubmanip},'dimCorrelation'} = [rho,p,h];
        
           %% TABLE
        % =================================
        subTable = table;
        fieldNames = {'condition','behavior'};
        for iField = 1:numel(fieldNames)
            varNames =  fieldnames(subdata.( fieldNames{iField}))';
            for iVar = 1:numel(varNames)
               subTable.(varNames{iVar}) = subdata.(fieldNames{iField}).(varNames{iVar})';
            end
        end
        dimTable = table( nominal(repmat( {submanipFieldName{iSubmanip}} ,height(subTable),1 )) ,'VariableNames',{'dimension'});
        subTable = [ dimTable  , subTable ];
        
        if isfield(data,'table')
            data.table = [ data.table ; subTable];
        else
            data.table = [ subTable];
        end
         
    end % file existence testing
    
        
end % submanip iteration
    
end
