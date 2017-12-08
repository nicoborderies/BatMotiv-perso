%% setGroups
%
% this script specify the groups of subjects to be analyzed
% as follows: 
%   - GROUP_XX:  a vector of subject id number
%   - groups:    a vector compiling all GROUP_XX
%   - groupNames: a cell list of string specifying the names of the
%                   different groups
%
% These variables are then used by the "analyze_batmotiv_group" function to
% look for the corresponding subject directory if organized as follows:
%   directory = 'blablabla\Bat-données\GROUP_XX\subid'
%
% NB. Customize this script to your own implementation, do not publish a public
%     modification
%
% Nicolas Borderies
% 11/2015


% control groups
    schizoControl = [ 13 , 16 , 17 , 19 , 20 , 25 , 31 , 32 , 33 , 37 , 40 , 45 , 46 , 47 , 48 , 53 ];
    youngControl = ([[101:106],[6001:6012]]);
    oldControl = ([7001:7007 , 7010:7022]);
    oldControl = ([7005:7007 , 7010:7026]); % control group for FTD study

    CONTROL = ([ schizoControl , youngControl , oldControl ]);
%     CONTROL = ([  oldControl ]);

    CITALOPRAM = ([1,3:12,14,16:23,25:28]);
    
    OPIOID = ([2 3 5 6 7 9 10:15 18 21:24 26 27 31 32 35:37]);
    OVERLOAD = ([ 1001 , 1:9 13:19 25:32 34 35 50:62]);
    
    IRM = [[8001:8015],[9000:9006,9008:9024]];
%     IRM = [[8001:8015],[9001:9004]];

% patient groups
    bvFTD = ([ 45 , 2005 ,2006, 2009 , 2010 ,2011 ,2012 ,2013 ,2014 ,2015 , 2016, 2017 , 2018 , 2020 , 2021 , 2022, 2023, 2024 , 2025, 2026 ]);
    % 2008 ??
    ACOM = ([ 1 2 3 4 6 7 8 10 12 13 15 21 23 25 27 31 32 33 34 35 37:42 44 52:58 60 61 62 64 65]);
    
    STEINERT = ([5,7,8,9,10,11,13,14,15,16,17,18,19,20]);
    SCHIZO = ([1:10 , 35 , 36 , 42 , 49 , 50 , 51 , 52 , 56 , 1001:1015 , 2001:2008,2010,3001:3007 ]); % 3022,3024:3028

    DEP = ([11,12,14,15,18,21,22,23,24,26,27,28,29,30,34,38,39,41,43,44,54,55]);

% compile
%     groups =  [ OVERLOAD ];
%     groupNames = {'OVERLOAD'};
%     groups =  [ CONTROL , bvFTD ];
%     groupNames = {'CONTROL','bvFTD'};