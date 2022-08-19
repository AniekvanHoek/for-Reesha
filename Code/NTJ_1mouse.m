%% Script plotting full hour session of individual mice, colored by behavioral events
% Aniek van Hoek
%
% Goal:     Plots full hour session in neural trajectory
%           Colors cue, water and ethanol indices
%           Also plots relative changes in direction for ethanol bouts
%
% Input:    Change mouse and date below      
%           
%
% Output:  Neural trajectories of full session, colored by events

clearvars; clc
close all;

cd = 'Reesha/MATLAB/cohort_9_calcium_imaging/Recent'; 
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
%% Mouse to analyse
date = "02082022";
mouse = "M5-2";

%% Loading neurons and indices of mouse
% Neural data: C & S  matrix from CNMFE analysis
name_neuron = sprintf('neuron/%s_%s_neuron.mat', date, mouse);
load(fullfile(cd, name_neuron)); 
allNeurons = neuron.C;

name_wateridx = sprintf('indices/%s_%s_water_idx.mat', date, mouse);
load(fullfile(cd,name_wateridx)); 

name_etohidx = sprintf('indices/%s_%s_ethanol_idx.mat', date, mouse);
load(fullfile(name_etohidx)); 

name_cueidx = sprintf('indices/%s_%s_cue_idx.mat', date, mouse);
load(fullfile(name_cueidx)); 

name_drinkingdata = sprintf('drinking_data/%s_%s_drinking_data.mat', date, mouse);
load(fullfile(cd,name_drinkingdata)); 

%% Change cue idx to every cue idx + 3 seconds
idx = find(drinking_data(2,:) == 1);
for i = 1:30
    drinking_data(2, idx(i):idx(i)+29) = 1;
end

%%  PCA
% Only mean center
pca_matrix = allNeurons;

% Mean center and scale PCA
for i = 1:size(allNeurons,1)
    mean_neuron = mean(allNeurons(i, :));
    std_neuron = std(allNeurons(i,:));
    
    pca_matrix(i,:) = (pca_matrix(i,:) - mean_neuron)/std_neuron;
end

[coef, ~, ~, ~, explained] = pca(pca_matrix(:,:)', 'Economy', false);
num_pcs = find(cumsum(explained)>=90, 1); 

data(:, :) = coef(:, 1: num_pcs)'*allNeurons(:,:);

 
%% Make relative plot (only etoh)
% Relative plot shows the direction in which the neural population response
% moves after a bout starts

for i = 1:size(ethanol_idx,2)
    startpoint = ethanol_idx(i);
    options = find(drinking_data(4,:) == 0);
    endpoint = options(find(options > ethanol_idx(i),1))-1;
    
    rel_data{i} = data(1:3,startpoint:endpoint) - data(1:3,startpoint-1);
end


smooth_window = [100 0];
d_marker_size = 20;

figure
for i = 1:size(rel_data,2)
    d1 = rel_data{i};
    d1 = smoothdata(d1, 2, 'gaussian', smooth_window);
   
    start = scatter3(d1(1, 1), d1(2, 1), d1(3, 1), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    finish = scatter3(d1(1, end), d1(2, end), d1(3, end), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
         
    start2 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 0]/255, '>', 'filled', 'Displayname', 'Start bout');
	finish2 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'Displayname', 'End bout');
  
    plot3(d1(1,:), d1(2,:), d1(3,:),'color',rand(1,3))
    hold on
end

    xlabel(sprintf("PC1: %0.2f%%", round(explained(1),2)))
    ylabel(sprintf("PC2: %0.2f%%", round(explained(2),2)))
    zlabel(sprintf("PC3: %0.2f%%", round(explained(3),2)))
    title(sprintf("Relative ethanol bouts of %s - %s", mouse, date));
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    

%% Neural trajectory plotting (full trial, line made invisible)
% Parameters
smooth_window = [10 0];
event_marker_size = 10;
line_opacity = 0;
line_width = 3;
errorbar_opacity = 1;
point_size = 5;
point_freq = 100;

scatter_plot_colors = {[0, 0, 0]/255}; 
line_plot_colors = {[211, 211, 211]/255};

N = 30; % Plot 30 differnt figures 
time = floor(size(data,2)/N);
counter = 0;
    s = get(0, 'ScreenSize');
    figure('Position', [0 0 s(3)  s(4)]);
    
for i = 1:N
    plot_data = {data(1:3,counter+1:counter+time)};

    w = waitforbuttonpress;
    plot1traj_3Dgaussian_AH_v2(plot_data, drinking_data(:,counter+1:counter+time), smooth_window, mouse, event_marker_size, ... %"choice none",
        {"Start", "Finish", "Ethanol", "Water", "Cue"}, line_plot_colors, line_opacity, line_width, scatter_plot_colors, point_size, point_freq)
    xlabel(sprintf("PC1: %0.2f%%", round(explained(1),2)))
    ylabel(sprintf("PC2: %0.2f%%", round(explained(2),2)))
    zlabel(sprintf("PC3: %0.2f%%", round(explained(3),2)))
    title(sprintf("Full neural trajectory of %s - %s, time %d", mouse, date, i));
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    counter = counter+time;
    if i == 1
        legend();
    end
end

