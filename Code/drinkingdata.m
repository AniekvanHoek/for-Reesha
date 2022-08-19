%% Getting drinking data of mice per session
% Aniek van Hoek
%
% Goal:     - This script outputs a matrix of dimensions 4 x # timepoints
%               with binary information on whether animal is drinking

% Input:    Change parameters below         
%
% Output:   Drinking data matrix
%           row 1:  time
%           row 2:  cue_idx
%           row 3:  water licks
%           row 4:  ethanol licks

clearvars; close all;
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
%% Parameters
date = '01082022';
cue = true; % true if cued 2BC experiments

%% Get drinking data 
for j = [4,5] % for both cages
    if j == 5
        for i = 1:3
            mousecage = j
            mousenum = i
            drinking_data = retrieving_drinking_data(date, mousecage, mousenum,cue);
            filename = sprintf("Reesha/MATLAB/cohort_9_calcium_imaging/drinking_data/%s_M%d-%d_drinking_data.mat", date, mousecage, mousenum);
            save(filename, 'drinking_data')
        end
    elseif j == 4
        for i = 2:3 % Skipping mouse 1
            mousecage = j
            mousenum = i
            drinking_data = retrieving_drinking_data(date, mousecage, mousenum,cue);
            filename = sprintf("Reesha/MATLAB/cohort_9_calcium_imaging/drinking_data/%s_M%d-%d_drinking_data.mat", date, mousecage, mousenum);
            save(filename, 'drinking_data')
        end

    end
end
    
