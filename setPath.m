function [root,vbadir,analysisdir,datadir,resultdir] = setPath
%
% this function specify the path required for the  "analyze_batmotiv_group" function
%
% NB. Customize this function to your own implementation, do not publish a public
%     modification
%
% Nicolas Borderies
% 11/2015


% define path
    % External Hard Drive 
%     hardDrive = 'F:'; % H
    hardDrive = 'H:'; % H
    icmDisk = 'B:';
    disk = icmDisk;
    
    root = [ disk '\nicolas.borderies\Projets scientifiques'];
    vbadir = [ disk '\nicolas.borderies\MATLAB\GitHub\VBA-toolbox'];
    
    batmotivdir = [root filesep 'batmotiv'];
    datadir=[batmotivdir filesep 'Bat-données'];
    analysisdir = [batmotivdir filesep 'Bat-analyse'];
    resultdir = [batmotivdir filesep 'Bat-résultats'];
    
% set path
%     addpath(genpath(datadir));
%     addpath(genpath(analysisdir));
%     addpath(genpath(resultdir));
%     addpath(genpath(vbadir));