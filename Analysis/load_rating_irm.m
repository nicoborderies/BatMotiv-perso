
function [data] = load_rating_irm(subjectDir,option)
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
        subManipList = {'ratingR_s','ratingR_im_s','ratingP_s', 'ratingE_s',...
                        'rating_money*gain','rating_money*loss','rating_grip'}; % to use as MANIP_NAME
        typeItemList={'rewardlist','rewardlist','punishlist', 'effortlist','stimlist','stimlist','stimlist'};            % to get list of item
        submanipFieldName={'writtenReward','visualReward', 'writtenPunishment', 'writtenEffort',...
                            'monetaryReward','monetaryPunishment','gripEffort'}; % to use as field name in data
        submanipAcronym = {'Rv','Ri','Pv','Ev','Rm','Pm','Eg'};
        nSubManip=length(submanipFieldName);
        data.misc = table(zeros(nSubManip,1),'RowNames',submanipFieldName,'VariableNames',{'nSession'});

    % session selection
        sessionList = [2:4]; 


%% ITERATION ACROSS TASKS
for iSubmanip = 1:nSubManip
    %% FILE EXTRACTION
    % ===========================================================================
    MANIP_NAME = subManipList{iSubmanip}; % result files should include manip 'manipName'
    
    fileList=dir([subjectDir filesep '*' MANIP_NAME '*.mat']);
    if isempty(fileList)
        warning(['No file found for ' MANIP_NAME 'in directory: ' subjectDir])
        isFileFound =0;
    else
        isFileFound =1;
    end
    
    
    
    if isFileFound == 1 % check file existence
        [fileList, nSession] = tools.load_file(subjectDir, MANIP_NAME, sessionList);
        data.misc{submanipFieldName{iSubmanip},'nSession'} = nSession;
        data.misc{submanipFieldName{iSubmanip},'acronym'} = {submanipAcronym{iSubmanip}};
        
        
        % concatenate sessions
        allData = vertcat(fileList.data);
        
        %% CONDITIONS
        % ===========================================================================
        subdata.condition.sessionNumber=[];
        for iSession = 1:nSession
            [nTrial,~]=size(fileList(iSession).data);
            subdata.condition.sessionNumber =   [subdata.condition.sessionNumber iSession * ones(1,nTrial)] ;
        end
        
        
        subdata.condition.trialNumber     = allData(:, 1)';
        permutation = horzcat(fileList.Perm) ;
        permutation = permutation(1:numel(subdata.condition.trialNumber));
        subdata.condition.itemNumber      = permutation  ;                                      % item number
        subdata.condition.itemSubtype     = cell(1,numel(subdata.condition.itemNumber ));
        
% subcatagorization of item
%         switch typeItemList{iSubmanip}
%             case 'rewardlist'
%                 subdata.condition.itemSubtype(subdata.condition.itemNumber>12)     =  repmat({'nonAlimentary'},1,sum(allData(:, 2)'>12)); 
%                 subdata.condition.itemSubtype(subdata.condition.itemNumber<=12)     = repmat({'alimentary'},1,sum(allData(:, 2)'<=12)); 
%             case 'punishlist' % be cautious with subtypes of punishments
%                 subdata.condition.itemSubtype(subdata.condition.itemNumber>13)     =  repmat({'nonSensory'},1,sum(allData(:, 2)'>13));
%                 subdata.condition.itemSubtype(subdata.condition.itemNumber<=13)     = repmat({'sensory'},1,sum(allData(:, 2)'<=13));
%             case 'effortlist'
%                 subdata.condition.itemSubtype(subdata.condition.itemNumber>12)     =  repmat({'cognitive'},1,sum(allData(:, 2)'>12));
%                 subdata.condition.itemSubtype(subdata.condition.itemNumber<=12)     = repmat({'motor'},1,sum(allData(:, 2)'<=12));
%             case 'stimlist'
%                 switch MANIP_NAME
%                     case 'rating_grip'
%                         subdata.condition.itemSubtype  =  repmat({'grip'},1,numel(allData(:, 2)));
%                     case {'rating_money*gain','rating_money*loss'}
%                         subdata.condition.itemSubtype  =  repmat({'money'},1,numel(allData(:, 2)));
%                 end
%                 
%         end
%         subdata.condition.itemSubtype = nominal(subdata.condition.itemSubtype);

    % preprocessing of item strings
    try
        itemList = horzcat(fileList(2:end).(typeItemList{iSubmanip}))  ;
        [l,~] = size(itemList);
        if l~=24; itemList = itemList';end
        
%         subdata.condition.itemName = nominal(itemList(subdata.condition.itemNumber)');      
        subdata.condition.itemName = nominal(itemList');        

%         countChar = @(s) numel(s);
%         countWord = @(s) numel(strfind(char(s),' '))+1;
        
        charNumberList = cell2mat(cellfun(@countChar,itemList,'UniformOutput',0));
        wordNumberList = cell2mat(cellfun(@countWord,itemList,'UniformOutput',0));

%         subdata.condition.itemCharNumber = charNumberList(subdata.condition.itemNumber)';
%         subdata.condition.itemWordNumber = wordNumberList(subdata.condition.itemNumber)';
        subdata.condition.itemCharNumber = charNumberList';
        subdata.condition.itemWordNumber = wordNumberList';
    catch
        subdata.condition.itemName = nominal(nan(size(subdata.condition.trialNumber)));   
        subdata.condition.itemCharNumber = nan(size(subdata.condition.trialNumber));   
        subdata.condition.itemWordNumber = nan(size(subdata.condition.trialNumber));   
    end
        
        
        % within-subject factor
        if ~isempty(load_withinSubFactor)
            [subdata.condition] =  load_withinSubFactor(subdata.condition,subjectDir,option);
        end
        
        %% BEHAVIOUR
        % ===========================================================================
        
        subdata.behavior.rating                 = allData(:, 2)';                  %
        subdata.behavior.ratingZValue           = zscore(subdata.behavior.rating);
        nR = numel(subdata.behavior.rating);
        subdata.behavior.ratingRank             = nan(1,nR);
        [i,j] = sort(subdata.behavior.rating);
        subdata.behavior.ratingRank(j) = [1:nR];
        
        
        subdata.behavior.reactionTime           = allData(:, 5)';       
        subdata.behavior.correctedReactionTime  = subdata.behavior.reactionTime./subdata.condition.itemCharNumber;       
        subdata.behavior.normalizedReactionTime = zscore(subdata.behavior.correctedReactionTime);
        
%         if isfield(fileList,'cursor')
%             subdata.behavior.cursor = fileList.cursor;
%             nC = numel(subdata.behavior.cursor);
%             pressNumber = zeros(1,nC);
%             directionNumber = zeros(1,nC);
%             pressSpeed = zeros(1,nC);
%             for iC = 1:nC
%                 c = subdata.behavior.cursor{iC};
%                 [pressNumber(iC),directionNumber(iC),pressSpeed(iC)] = count_pressnumber(c);
%             end
%             subdata.behavior.pressNumber = pressNumber;
%             subdata.behavior.directionNumber = directionNumber;
%             subdata.behavior.pressSpeed = pressSpeed;
%         end
        
        %% MISC
        % ===========================================================================
        data.misc{submanipFieldName{iSubmanip},'fileName'} = {fileList.fileName};
        try
        data.misc{submanipFieldName{iSubmanip},'itemList'} = {fileList.([typeItemList{iSubmanip}])};
        end

        
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
        
    end
end


end

function c = countChar(s)
    c = numel(s);
end

function c = countWord(s)
    c = numel(strfind(char(s),' '))+1;
end