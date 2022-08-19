function [water_bouts, neural_data_water, etoh_bouts, neural_data_etoh] = get_epochs(drinking_data, neural_activity, min_distance)
%  Aniek van Hoek
%
% Goal:     Get bouts of water and substance that are far enough apart

% Input:    
%   drinking_data:      Matrix with binary information on drinking
%   neural_activity:    Neural activity matrix per trial
%   min_distance:       how far bouts have to be apart
%
% Output:
%   water_bouts:        Cell array of water bouts timeframe indices
%   neural_data_water:  Neural data corresponding to timeframes bouts
%   etoh_bouts:         Cell array of alcohol bouts timeframe indices
%   neural_data_etoh:   Neural data corresponding to timeframes bouts
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
%% Water
water = drinking_data(3,:);
wtr_idx = find(water);

% Make cell array of different bouts indices for water    
idx = find(diff(wtr_idx)~= 1) +1;
for i = 1:size(idx,2)

   if i ~= size(idx,2)
        water_bouts{i+1} = wtr_idx(idx(i)):wtr_idx(idx(i+1)-1); % For last bout
   else 
        water_bouts{i+1} = wtr_idx(idx(i)):wtr_idx(end); % All other bouts
   end
   
end
water_bouts{1} = wtr_idx(1):wtr_idx(idx(1)-1);


% Delete bouts that are too close to each other
previous = 1;
for i = 1:size(water_bouts,2)-1
    try
        if water_bouts{i+1}(1) - min_distance <= water_bouts{i}(end) 
            water_bouts{i+1}= [];
        end
        
    catch ME
        if water_bouts{i+1}(1) - min_distance <= water_bouts{i-previous}(end) 
            water_bouts{i+1}= [];
            previous = previous+1;
        end

    end
end  

water_bouts = water_bouts(~cellfun('isempty',water_bouts));
    
% Get neural epochs
for m = 1:size(water_bouts,2)
    neural_data_water{m} = neural_activity(:,water_bouts{m});
end


%% Ethanol
ethanol = drinking_data(4,:);
etoh_idx = find(ethanol);
% Make cell array of different bouts indices for water    
idx = find(diff(etoh_idx)~= 1) +1;
for i = 1:size(idx,2)

   if i ~= size(idx,2)
        etoh_bouts{i+1} = etoh_idx(idx(i)):etoh_idx(idx(i+1)-1);
   else 
        etoh_bouts{i+1} = etoh_idx(idx(i)):etoh_idx(end);
   end
   
end
etoh_bouts{1} = etoh_idx(1):etoh_idx(idx(1)-1);


% Delete bouts that are too close to each other
previous = 1;
for i = 1:size(etoh_bouts,2)-1
    try
        if etoh_bouts{i+1}(1) - min_distance <= etoh_bouts{i}(end) 
            etoh_bouts{i+1}= [];
        end
        
    catch ME
        if etoh_bouts{i+1}(1) - min_distance <= etoh_bouts{i-previous}(end) 
            etoh_bouts{i+1}= [];
            previous = previous+1;
        end

    end
end  

etoh_bouts = etoh_bouts(~cellfun('isempty',etoh_bouts));
    
% Get neural epochs
for m = 1:size(etoh_bouts,2)
    neural_data_etoh{m} = neural_activity(:,etoh_bouts{m});
end

end
