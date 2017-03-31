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
    
    root = [ disk '\nicolas.borderies\projets'];
    vbadir = [ disk '\nicolas.borderies\MATLAB\GitHub\VBA-toolbox'];
    
    batmotivdir = [root filesep 'batmotiv'];
    datadir=[batmotivdir filesep 'données'];
    analysisdir = [batmotivdir filesep 'code.public'];
    resultdir = [batmotivdir filesep 'resultats'];
    
% set path
%     addpath(genpath(datadir));
%     addpath(genpath(analysisdir));
%     addpath(genpath(resultdir));
%     addpath(genpath(vbadir));