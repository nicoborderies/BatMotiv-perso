
function [data] = load_choice_irm(subjectDir,option,resultStructure)
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
    subManipList={'choiceR_s','choiceR_im_s','choiceP_s', 'choiceE_s'}; % to use as MANIP_NAME
    typeItemList={'rewardlist','rewardlist','punishlist', 'effortlist'};            % to get list of item
    submanipAcronym = {'Rv','Ri','Pv','Ev'};
     itemCategories = {'alimentary','nonAlimentary';
                       'alimentary','nonAlimentary';
                       'nonSensory','sensory';
                       'cognitive','motor'};

    submanipFieldName={'writtenReward','visualReward', 'writtenPunishment', 'writtenEffort'}; % to use as field name in data
    nSubManip=length(submanipFieldName);
    data.misc = table(zeros(nSubManip,1),'RowNames',submanipFieldName,'VariableNames',{'nSession'});

    % session selection
        sessionList = [2:4]; 



%% ITERATION ACROSS TASKS
for iSubmanip = 1:nSubManip
    %% FILE EXTRACTION
    % ===========================================================================
    MANIP_NAME = subManipList{iSubmanip};
    
    
    fileList=dir([subjectDir filesep '*' MANIP_NAME '*.mat']);
    if isempty(fileList)
        warning(['No file found for ' MANIP_NAME ' in directory ' subjectDir])
        isFileFound =0;
    else
        isFileFound =1;
    end
    
%     sessionList = {fileList(2:4).name};
    
    if isFileFound == 1
       [fileList, nSession] = tools.load_file(subjectDir, MANIP_NAME, sessionList);
        data.misc{submanipFieldName{iSubmanip},'nSession'} = nSession;
        data.misc{submanipFieldName{iSubmanip},'acronym'} = {submanipAcronym{iSubmanip}};
        
                %% MISC
        % ===========================================================================
         % concatenate sessions
        allData = vertcat(fileList.data);
        
         % --- correction for erroneous saving of rating*Dim*summary --- %
%         get_correctRating_choice;
        data.misc{submanipFieldName{iSubmanip},'fileName'} = {fileList.fileName};

         % rating summaries
%          switch submanipFieldName{iSubmanip}
%             case {'writtenReward'}
%                   data.misc{submanipFieldName{iSubmanip},'ratingSummary'} = {fileList.ratingRsummary};
%             case {'visualReward'}
%                   data.misc{submanipFieldName{iSubmanip},'ratingSummary'} = {fileList.ratingRsummary};
%             case {'writtenPunishment'}
%                   data.misc{submanipFieldName{iSubmanip},'ratingSummary'} = {fileList.ratingPsummary};
%             case {'writtenEffort'}
%                   data.misc{submanipFieldName{iSubmanip},'ratingSummary'} = {fileList.ratingEsummary};
%       
%          end
         
        
        %% CONDITIONS
        % ===========================================================================
        subdata.condition.sessionNumber=[];
        for iSession = 1:nSession
            [nTrial,~]=size(fileList(iSession).data);
            subdata.condition.sessionNumber =   [subdata.condition.sessionNumber iSession * ones(1,nTrial)] ;
        end
        
         
        signDeltaValue =  strcmp(submanipFieldName{iSubmanip},'writtenReward')...
                        + strcmp(submanipFieldName{iSubmanip},'visualReward')...
                        - strcmp(submanipFieldName{iSubmanip},'writtenPunishment')...
                        - strcmp(submanipFieldName{iSubmanip},'writtenEffort'); 
                    
        subdata.condition.trialNumber          = allData(:, 1)';
        subdata.condition.leftItemNumber       = allData(:, 4)';                         % left item number (index, see rating manip to get actual items)
        subdata.condition.rightItemNumber      = allData(:, 6)';                         % right item number (index, see rating manip to get actual items)
        subdata.condition.leftRatingValue      = allData(:, 5)';                         % left item value (from rating)
        subdata.condition.rightRatingValue     = allData(:, 7)';                         % right item value (from rating)
        subdata.condition.differenceRatingValue      = (allData(:, 7)' - allData(:, 5)'); % difference rating (right - left)
%         itemsCategory = {'leftItemCategory','rightItemCategory'};
%         itemsNumber = {'leftItemNumber','rightItemNumber'};
%         itemsName = {'leftItemName','rightItemName'};
%         preferedCategory = {'leftPreferedCategory','rightPreferedCategory'};
%         
%         [itemList] = get_itemList_choice(subjectDir,submanipFieldName{iSubmanip});
%         [l,~] = size(itemList);
%         if l~=24; itemList = itemList';end      
%         for iT=1:2
%             subdata.condition.(itemsCategory{iT}) = convert2category(subdata.condition.(itemsNumber{iT}),typeItemList{iSubmanip});
%             subdata.condition.(itemsCategory{iT}) = nominal(subdata.condition.(itemsCategory{iT}));
%             subdata.condition.(itemsName{iT}) =  nominal(itemList(subdata.condition.(itemsNumber{iT}))');
%         end
%         meanCat = tools.tapply(subdata.condition.leftRatingValue, {subdata.condition.leftItemCategory}, @nanmean);  
%         [m,index] = max(meanCat);
%         pref =  nominal(itemCategories{iSubmanip,index});
%         for iT=1:2
%             subdata.condition.(preferedCategory{iT}) = ismember(subdata.condition.(itemsCategory{iT}),pref);
%         end
        
        % zscore ratings
%         mu=mean(data.misc{submanipFieldName{iSubmanip},'ratingSummary'}{:});
%         sigma=std(data.misc{submanipFieldName{iSubmanip},'ratingSummary'}{:});
%         v_min  = min(data.misc{submanipFieldName{iSubmanip},'ratingSummary'}{:}); 
%         v_range = range(data.misc{submanipFieldName{iSubmanip},'ratingSummary'}{:}); 
%         subdata.condition.rightZValue = (subdata.condition.rightRatingValue-mu)/sigma;
%         subdata.condition.leftZValue  = (subdata.condition.leftRatingValue-mu)/sigma;
%         subdata.condition.differenceZValue = [subdata.condition.rightZValue -  subdata.condition.leftZValue];
%         
%         subdata.condition.rightRangeValue = (subdata.condition.rightRatingValue - v_min)/v_range;
%         subdata.condition.leftRangeValue  = (subdata.condition.leftRatingValue - v_min)/v_range;
        
        % consensual ratings
%         try
%             cd ..
%             cd ..
%             load('consensualValue.mat');
%             cd(subjectDir);
%             [~,ind] = ismember(subdata.condition.leftItemName',consensus.itemName);
%             subdata.condition.leftConsensusValue = consensus.value(ind)';
%             [~,ind] = ismember(subdata.condition.rightItemName',consensus.itemName);
%             subdata.condition.rightConsensusValue = consensus.value(ind)';
% 
%             % search sub index
%                 groupName = option.design.group;
%                 subind = strfind(subjectDir,'sub');
%                 subid = str2num((subjectDir(subind+3:end)));
%                 ind = find(scaling.group==groupName &...
%                            scaling.subject==subid &...
%                            scaling.dimension==submanipFieldName{iSubmanip})  ;
%   
%                 gain = scaling.gain(ind);
%                 offset = scaling.offset(ind);
%                 
%             subdata.condition.leftPersonnalValue = (subdata.condition.leftRatingValue./100)./gain - offset; 
%             subdata.condition.rightPersonnalValue = (subdata.condition.rightRatingValue./100)./gain - offset; 
%         catch
%             subdata.condition.leftPersonnalValue = nan(1,numel(subdata.condition.leftRatingValue)); 
%             subdata.condition.rightPersonnalValue = nan(1,numel(subdata.condition.leftRatingValue)); 
%         end
        
         % within-subject factor
        if ~isempty(load_withinSubFactor)
            [subdata.condition] =  load_withinSubFactor(subdata.condition,subjectDir,option);
        end
        
        %% BEHAVIOUR
        % ===========================================================================
        
        subdata.behavior.RT  = allData(:, 8)';                          % RT
        subdata.behavior.sideChoice  = allData(:, 9)';                  % choice side: [1] right [-1] left
%         subdata.behavior.preferedChoice(subdata.behavior.sideChoice==-1)  = isequal(subdata.condition.leftItemCategory(subdata.behavior.sideChoice==-1) , subdata.condition.leftPreferedCategory(subdata.behavior.sideChoice==-1)); 
%         subdata.behavior.preferedChoice(subdata.behavior.sideChoice==1)  = isequal(subdata.condition.rightItemCategory(subdata.behavior.sideChoice==-1),subdata.condition.rightPreferedCategory(subdata.behavior.sideChoice==-1)); 

  
        
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



function cat = convert2category(number,dimension)
    cat = cell(1,numel(number));
        switch dimension
            case 'rewardlist'
                cat(find(number>12)) = repmat({'nonAlimentary'},1,numel(find(number>12))); 
                cat(find(number<=12)) = repmat({'alimentary'},1,numel(find(number<=12))); 
            case 'punishlist'
                cat(find(number>13)) = repmat({'nonSensory'},1,numel(find(number>13))); 
                cat(find(number<=13)) = repmat({'sensory'},1,numel(find(number<=13))); 
            case 'effortlist'
                cat(find(number>12)) = repmat({'cognitive'},1,numel(find(number>12))); 
                cat(find(number<=12)) = repmat({'motor'},1,numel(find(number<=12))); 

        end
end





