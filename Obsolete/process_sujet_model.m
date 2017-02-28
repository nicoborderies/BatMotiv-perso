function result = process_sujet_model(subjectData,modelName, option)
% PRECESS_SUJET_MODEL process subject's data with the modelName processor
% result = process_sujet_model(subjectData,modelName, optionList)
% IN
%      - subjectData : subject's structure
%      - modelName   : name of the processor
%      - option      : option to pass to the processor
% OUT
%      - result      


% default arg
if nargin < 3, option = []; end

% check if manip loader exists
processorName = ['+processors/process_data_' modelName];
if exist(processorName,'file') ~= 2
        error('*** Cannot find the function %s, check the manipName.',processorName);
end

% load datafile
% double-loop allowing to load data even if some fields are missing !
processor_fun = str2func(['processors.process_data_' modelName]);
result=feval(processor_fun,subjectData,option);

end

