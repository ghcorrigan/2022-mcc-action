%% Midcingulate cortex during stopping
% 2023-mcc-action, S P Errington, August 2023

%% Dev notes:
%   2023-08-26: Continuing development.

%% Setup workspace
% Clear workspace
clear all; close all; clc; beep off; warning off;

% Define paths & key directories
system_id = 'robarts';
dirs = get_dirs_mcc(system_id);
getColors_mcc; % Define colorschemes

%% Data extraction
% Get recording log
ephysLog = importOnlineEphysLogMaster; 

% Get dataset map
dajo_datamap = load_datamap(fullfile(dirs.root,'data')); 


%% Analysis scripts

midcingulate_primary

dev_dmfc_analysis


