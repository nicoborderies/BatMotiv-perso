%% analyze_batmotiv_opioid_2level_gripox

% reset
clc;
clear all;
close all;

%% Prepare data analysis
%-------------------------------

% load dataset
%%% design table
datadir = 'B:\nicolas.borderies\projets\batmotiv\données\OPIOID';
cd(datadir);
f = load('subject_table_opioid.mat');
design = f.tab;  
% processed dataset
datadir = 'B:\nicolas.borderies\projets\batmotiv\resultats\OPIOID';
codedir = 'B:\nicolas.borderies\projets\batmotiv\code.perso';
cd(datadir);
analysisname = ['batmotiv_analysis_OPIOID_gripox_all_30_5_2017.mat'];
load(analysisname);
cd(codedir);


        
%% Data analysis
%-------------------------------

% define analysis parameters
%%% reformat
data = [groupData.OPIOID.subject ];
result = [groupResult.OPIOID.subject ];
%%% lists
dimensionList = {'writtenReward','visualReward','writtenPunishment','writtenEffort',...
                            'monetaryReward','monetaryPunishment','gripEffort',...
                            'RewardEffort','PunishmentEffort','RewardPunishment',...
                            'Gain','Loss','GainLoss','GainEmotion'};
taskList =  {'rating','choice','weight','discount',...
                    'grip','gripIAPS','gripAccu','mental','learning'};
treatmentList = {'naloxone','placebo','morphine'};
orderGrip = design.order_taskGrip(design.selection==1);
orderRating = design.order_rating(design.selection==1);
orderWeight = design.order_weight(design.selection==1);
weight = design.weight(design.selection==1);

ntrt = numel(treatmentList);
nsub = numel(data);
%%% display
col = { [0 0 1]*0.75 , [1 1 1]*0.5 , [1 0 0]*0.75 };

%% compare oximetric statistics between treatments


hr = nan(ntrt,nsub);
hrv = nan(ntrt,nsub);


for isub=1:nsub % subject loop
    % test completeness of oximetric data
   if isfield(result{isub}.grip,'oximeter') 
       if numel(result{isub}.grip.oximeter)==3
          
           % session-to-treatment mapping
           datatab = data{isub}.grip.table  ; 
           trt = datatab.treatment;
           trt = removecats(trt,'0');
           trt = reordercats(trt,treatmentList);
           sess = datatab.sessionNumber;
           trt_by_sess = splitapply(@unique,trt,sess);

           % extract
           for isess=1:3 % session loop
               stat = result{isub}.grip.oximeter(isess).oxistat;
               oxidata = result{isub}.grip.oximeter(isess).oxi;
               if ~isempty(stat)
    %                hr(double(trt_by_sess(isess)),isub) = stat.mean_hr;
                   hr(double(trt_by_sess(isess)),isub) = median(oxidata.hr);
                   hrv(double(trt_by_sess(isess)),isub) = stat.hrv;
               end                
           end
       end
   end
end



