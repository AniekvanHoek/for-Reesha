% Aniek van Hoek
%
% Goal:     Script saves C-matrix for Felix

clc; clearvars;
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));

%%% Retrieving of Data %%%
date = ["01082022", "01232022", "02082022", "02112022", "02282022", "03012022"];
mouse_full_load = ["M4-2", "M4-3", "M5-1","M5-2","M5-3"];
mouse_full_save = ["4-2", "4-3", "5-1","5-2","5-3"];

for j = 1:6
for x = 1:5
mouse_load = mouse_full_load(x);
mouse_save = mouse_full_save(x);
%% Getting Neural data
name_neuron = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/neuron/%s_%s_neuron.mat', date(j), mouse_load);
load(name_neuron); 
C_matrix = neuron.C;


name_save = sprintf('Reesha/MATLAB/cohort_9_calcium_imaging/Recent/C_matrices/%s_%s_neuron', date(j), mouse_load);
save(name_save, 'C_matrix');
end
end