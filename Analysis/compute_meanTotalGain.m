%% compute_meanTotalGain
% compute average total gain won by subject from a control group
%
% Nicolas Borderies- March 2017


% prepare variables
gripGain = [];
stroopGain = [];
learningGain = [];

% load
groupdir = 'B:\nicolas.borderies\projets\batmotiv\données\CONTROL';
files = dir(groupdir);


for i=1:numel(files)
    
    fprintf('loading file %d /%d \n',i,numel(files));
    
    % detect subject dir
    isSub = ~isempty(strfind(files(i).name,'sub'));
    if isSub 
       subdir = [groupdir filesep files(i).name ];
       cd(subdir);
       
       % extract grip gain
       file = dir([subdir filesep '*_gripRP*.mat']);
       if ~isempty(file) && numel(file)==1
          load(file.name);
%           gripGain = [gripGain; data(end,7)]; % total gain
          gripGain = [gripGain; mean(data(:,5)) ]; % performance
       end
       
      % extract stroop gain
       file = dir([subdir filesep '*mental*.mat']);
       if ~isempty(file) && numel(file)==1
          load(file.name);
%           stroopGain = [stroopGain; data(end,7)]; % total gain
          stroopGain = [stroopGain;  mean(data(:,5)) ]; % performance
       end
       
       % extract learning gain
       file = dir([subdir filesep '*learning*.mat']);
       if ~isempty(file) && numel(file)==3
          subdata = [];
          for is=1:3
            load(file(is).name);
            subdata = [subdata;data];
          end
          learningGain = [learningGain; sum(subdata(:,10))]; % total gain
%           learningGain = [learningGain; mean((subdata(subdata(:,3)~=2,8)+1)/2) ]; % performance
       end
       
    end
end

% display
disp(nanmean(gripGain));
disp(nanmean(stroopGain));
disp(nanmean(learningGain));