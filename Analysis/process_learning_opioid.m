function [result,data,option] = process_learning_opioid(data,option)
%% INITIALIZE 

       init_process_task;

    nBin = 6;
    
    
%% ANALYSES
% ===========================================================================
     % iteration across withinSubFactor
    factor = data.table.(withinSubFactor{1});
    levels = unique(factor); if ~iscell(levels);levels = num2cell(levels);end
        nLevel = numel(unique(factor));

    
        % pre-allocation
        aNum = [];
        for iLevels = 1:nLevel
         a = ones(nAnalysis,1)*iLevels ;
         aNum = [aNum ; a];
        end
        descriptive = table(  unique(factor) ,'VariableNames',withinSubFactor);
        inferential = table(  unique(factor) ,'VariableNames',withinSubFactor);
        misc 	    = table(  repmat([1:nAnalysis]',nLevel,1),aNum ,'VariableNames',{'analysisNumber',withinSubFactor{:}});
        model = struct; 
    
    for iLevels = 1:numel(levels)
        model(iLevels).(withinSubFactor{1}) = iLevels;
        
     
           % data selection (with respect to the task & factor)
                index = find(  factor == levels{iLevels} );
                design = data.table(index,:);
                
                % performance index
                y = tools.tapply(design.isOptimalChoice, {design.pairValence}, @nanmean);
                inferential{iLevels,['correct_Gain']} = y(3); 
                inferential{iLevels,['correct_Loss']} = y(1); 
                
                % repetition index
                y = tools.tapply(design.exploitationChoice, {design.pairValence}, @nanmean);
                inferential{iLevels,['repetition_Neutral']} = y(2);
                
                % WS/LS index
                y = tools.tapply(design.exploitationChoice, {design.previous_outcome,design.previous_bestOutcome}, @nanmean);
                inferential{iLevels,['repetition_previousLoss']} = y(1,1);
                inferential{iLevels,['repetition_previousNeutral_Loss']} = y(2,1);
                inferential{iLevels,['repetition_previousNeutral_Gain']} = y(2,2);
                inferential{iLevels,['repetition_previousGain']} = y(3,2);
                
                % learning index
                nt = quantileranks(design.trialNumberByValence,6);
                y = tools.tapply(design.isOptimalChoice, {design.pairValence,nt}, @nanmean);
                dy = nanmean(diff(y,1,2),2);
                inferential{iLevels,['dcorrect_Gain']} = dy(3); 
                inferential{iLevels,['dcorrect_Loss']} = dy(1); 
                
                confirmation = design.outcomePredictibility;
                valence = design.pairValence;
                block = design.blockNumber;
                nt = design.trialNumberByValence;
                confirmation_t = nan(size(confirmation));
                for iv = [-1 1]
                    for b = unique(block)'
                        confirmation_t(valence==iv & block==b) = [ NaN ; confirmation(valence==iv & block==b & nt~=max(nt)) ];
                    end
                end
                y = tools.tapply(design.isOptimalChoice, {design.pairValence,confirmation_t}, @nanmean);
                inferential{iLevels,['correct_infirmed_Loss']} = y(1,1);
                inferential{iLevels,['correct_confirmed_Loss']} = y(1,2);
                inferential{iLevels,['correct_infirmed_Gain']} = y(3,1);
                inferential{iLevels,['correct_confirmed_Gain']} = y(3,2);

                % optimality index
                glm = fitglm(design.previous_meanDiffOutcome,design.isOptimalChoice,...
                       'correct ~ -1 + dv',...
                       'VarNames',{'dv','correct'},...
                       'Link','logit','Distribution','Binomial');
                beta = glm.Coefficients.Estimate;
                inferential{iLevels,['beta_DV']} = beta; 
                
                % working memory index
                nt = design.trialNumber;
                valence = design.pairValence;
                interlapstime = nan(size(nt));
                for iv = [-1 1]
                    for b = unique(block)'
                        interlapstime(valence==iv & block==b) = [ NaN ; diff(nt(valence==iv & block==b)) ];
                    end
                end
                y = tools.tapply(design.isOptimalChoice, {interlapstime}, @nanmean);
                beta = glmfit(interlapstime,design.isOptimalChoice,'normal');
                inferential{iLevels,['beta_elapsedTime']} = beta(2);


                
                
                


     % model-based inverse inference (VBA)
      for iAnalysis = 1:nAnalysis
             
             fprintf('model:  %d / %d  \n',iAnalysis,nAnalysis);

             % set analysis
             [ ~,option ] = set_analysis( option,iAnalysis);

            % model
             [model(iAnalysis + (iLevels-1)*nAnalysis).vba,option] = invert_learning_opioid(design,option);
              vba = model(iAnalysis + (iLevels-1)*nAnalysis).vba;
              misc{iAnalysis + (iLevels-1)*nAnalysis, ['logE_learning']} = vba.logE;
      end
        selectedModel = option.analysis.learning.selectedModel ; % full model as a prior to estimate parameters
        
        % concatenate design matrix with model predictions
            predictions = model(selectedModel + (iLevels-1)*nAnalysis).vba.predictions;
            data.table.predicted_isOptimalChoice(index,1)  = predictions.isOptimalChoice';
            design = data.table(index,:);
            
        % extract parameters
        parameters = model(selectedModel + (iLevels-1)*nAnalysis).vba.parameters;
        paramNames = fieldnames(parameters);
        for iP = 1:numel(paramNames)
              inferential{iLevels,[ paramNames{iP} ]}  = [ parameters.(paramNames{iP})(1) ];
        end
        inferential{iLevels,['BCA']}  = [model(selectedModel + (iLevels-1)*nAnalysis).vba.BCA ];
      
     end

   
    
    % save result
        store_task_result;

 
   
end




