%% analyze_motiscan_opioid_psychometry

% reset
% ----------
clear all;
close all;

%% preparation
% folders definitions
datadir = 'B:\nicolas.borderies\projets\batmotiv\données\OPIOID';
codedir = 'B:\nicolas.borderies\projets\batmotiv\code.public\Bat-analyse';
subcodedir = 'B:\nicolas.borderies\projets\batmotiv\code.perso\Analysis\motiscan_opioid';
resultdir = 'B:\nicolas.borderies\projets\batmotiv\resultats\OPIOID';
addpath(genpath(codedir));
addpath(genpath(subcodedir));

% file definitions
resultname = 'motiscan-opioid-1level';
filenames = {'MBB_battery_psychometry','norris_score'};


% load experimental design data
cd(datadir);
design = load('subject_table_opioid.mat');

% subject definition
sub4analysis = design.tab.n_inclusion(find(design.tab.selection==1))';
nsub = numel(sub4analysis);
nsession = 3;

%% data aggregation

% data preparation
psychometry = cell(1,nsub);

for isub = sub4analysis % subject loop
    
    % data preparation
    psychometry{isub} = struct('session',{1,2,3});
    
    % display
    i = find(sub4analysis==isub);
    fprintf('analysis: subject %d / %d \n',i,nsub);
    
    % folder definition
    subdir = [datadir filesep 'sub' num2str(isub) ];
    cd(subdir);
    
    % load & format data
    % - 1. other scales
    f = dir([ filenames{1} '_sub' num2str(isub) '*.mat' ]);
    nfiles = numel(f);
    for ifile = 1:nfiles
       d = load(f(ifile).name) ;
       psychometry{isub}(ifile).BARATT = d.psychometry.BARATT ;
       psychometry{isub}(ifile).BIS = d.psychometry.BIS ;
       psychometry{isub}(ifile).LAY = d.psychometry.LAY ;
       psychometry{isub}(ifile).STARKSTEIN = d.psychometry.STARKSTEIN ;
    end   
    
    % - 2. norris scale 
    f = dir([ filenames{2} '_sub' num2str(isub) '*.mat' ]);
    nfiles = numel(f);
    for ifile = 1:nfiles
       d = load(f(ifile).name) ;
       psychometry{isub}(ifile).NORRIS = d ;
    end
    
       
end

% table creation
psychoTab = table(repelem(sub4analysis,1,nsession)','VariableNames',{'subid'});
psychoTab.session = nan(nsub*nsession,1);
psychoTab.treatment = categorical(repmat({''},nsub*nsession,1));

for isub = sub4analysis
    psychoTab.session(psychoTab.subid==isub) = [ 1 2 3]';
    psychoTab.treatment(psychoTab.subid==isub) = categorical([ design.tab.treatment_sess1(design.tab.n_inclusion==isub),...
                                                   design.tab.treatment_sess2(design.tab.n_inclusion==isub),...
                                                    design.tab.treatment_sess3(design.tab.n_inclusion==isub)]');
    for isess = 1:nsession
        if ~isempty(psychometry{isub}(isess))
            try
                psychoTab.apathy(psychoTab.subid==isub & psychoTab.session==isess) = psychometry{isub}(isess).STARKSTEIN.total;
                psychoTab.impulsivity(psychoTab.subid==isub & psychoTab.session==isess) = psychometry{isub}(isess).BARATT.total;
                psychoTab.activation(psychoTab.subid==isub & psychoTab.session==isess) = sum([psychometry{isub}(isess).BIS.normalized.BAS_drive ,...
                                                                                              psychometry{isub}(isess).BIS.normalized.BAS_funSeeking,...
                                                                                              psychometry{isub}(isess).BIS.normalized.BAS_rewardResponsiveness ]);
                psychoTab.inhibition(psychoTab.subid==isub & psychoTab.session==isess) = psychometry{isub}(isess).BIS.normalized.BIS;
                psychoTab.procrastination(psychoTab.subid==isub & psychoTab.session==isess) = psychometry{isub}(isess).LAY.total;
                psychoTab.vigilance(psychoTab.subid==isub & psychoTab.session==isess) = psychometry{isub}(isess).NORRIS.mean_score;
            catch
                psychoTab.apathy(psychoTab.subid==isub & psychoTab.session==isess) = NaN;
                psychoTab.impulsivity(psychoTab.subid==isub & psychoTab.session==isess) =  NaN;
                psychoTab.activation(psychoTab.subid==isub & psychoTab.session==isess) =  NaN;
                psychoTab.inhibition(psychoTab.subid==isub & psychoTab.session==isess) =  NaN;
                psychoTab.procrastination(psychoTab.subid==isub & psychoTab.session==isess) =  NaN;
                psychoTab.vigilance(psychoTab.subid==isub & psychoTab.session==isess) = NaN;
            end

        end
    end
end

psychoTab.treatment =  reordercats(psychoTab.treatment,{'naloxone','placebo','morphine'});

%% analysis

scale_name = {'apathy','impulsivity','activation','inhibition','procrastination','vigilance'};
for is = 1:numel(scale_name)
    
    scale = scale_name{is};
    disp(scale);
    g = findgroups(psychoTab.treatment);
    fname = @(x) {identity(x)};
    Y = splitapply(fname,psychoTab.(scale),double(psychoTab.treatment));
    Y = cell2mat(Y');
    Y = Y - nanmean(Y,2);
    [h,p,~,stat] = ttest(Y(:,3),Y(:,1));
    fprintf('opioid t-test: h = %d , p = %d , t = %d \n',h,p,stat.tstat);
    fprintf('\n');
    
end



