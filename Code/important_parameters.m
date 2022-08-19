%% Find important parameters in decoding results
% Aniek van Hoek
%
% Goal:     Scripts reads out the weights of neurons for the five best
%           decoding results (water vs substance) using co-registered 
%           neurons and a linear SVM.           
%
% Input:    Automated, choose sucrose or ethanol results
%
% Output:   Distribution plots of parameters
%           List of important neurons per session decoding

close all; 
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
sucrose = false;

%% Get data
if sucrose
    par_ses1 = readmatrix("Reesha/MATLAB/cohort_9_calcium_imaging/Recent/decoding/parameter_Sets.xlsx", 'Sheet','Pre-iso suc', 'Range','A1:E119');
    par_ses2 = readmatrix("Reesha/MATLAB/cohort_9_calcium_imaging/Recent/decoding/parameter_Sets.xlsx", 'Sheet','iso suc', 'Range','A1:E119');
    par_ses3 = readmatrix("Reesha/MATLAB/cohort_9_calcium_imaging/Recent/decoding/parameter_Sets.xlsx", 'Sheet','post-iso suc', 'Range','A1:E119');

else
    par_ses1 = readmatrix("Reesha/MATLAB/cohort_9_calcium_imaging/Recent/decoding/parameter_Sets.xlsx", 'Sheet','Pre-iso etoh', 'Range','A1:E131');
    par_ses2 = readmatrix("Reesha/MATLAB/cohort_9_calcium_imaging/Recent/decoding/parameter_Sets.xlsx", 'Sheet','Iso etoh', 'Range','A1:E131');
    par_ses3 = readmatrix("Reesha/MATLAB/cohort_9_calcium_imaging/Recent/decoding/parameter_Sets.xlsx", 'Sheet','post-iso etoh', 'Range','A1:E131');
end
    
%% Finding important neurons per session
% Session 1
% Distribution
figure
histogram(par_ses1);
title("session 1 etoh")

% Find boundaries of 2SD from mean
bounds_ses1 = [mean(par_ses1,'all') - 2*std(par_ses1,[],'all'),mean(par_ses1,'all')+ 2*std(par_ses1,[],'all')];

% Find indices of neurons who are in outside 95%
[row1,col1] = find(par_ses1>bounds_ses1(2));
[row2,col2] = find(par_ses1<bounds_ses1(1));

% Keep only unique neurons
final_ses1 = vertcat(unique(row1), unique(row2));



% Session 2
figure
histogram(par_ses2);
title("session 2 etoh")

bounds_ses2 = [mean(par_ses2,'all') - 2*std(par_ses2,[],'all'),mean(par_ses2,'all')+ 2*std(par_ses2,[],'all')];

[row1_ses2,col1_ses2] = find(par_ses2>bounds_ses2(2));
[row2_ses2,col2_ses2] = find(par_ses2<bounds_ses2(1));

final_ses2 = vertcat(unique(row1_ses2), unique(row2_ses2));



% Session 3
figure
histogram(par_ses3);
title("session 3 etoh")

bounds_ses3 = [mean(par_ses3,'all') - 2*std(par_ses3,[],'all'),mean(par_ses3,'all')+ 2*std(par_ses3,[],'all')];

[row1_ses3,col1_ses3] = find(par_ses3>bounds_ses3(2));
[row2_ses3,col2_ses3] = find(par_ses3<bounds_ses3(1));

final_ses3 = vertcat(unique(row1_ses3), unique(row2_ses3));