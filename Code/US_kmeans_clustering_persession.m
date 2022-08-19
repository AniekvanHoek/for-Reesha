%% Clustering of US responses using k-means
% Aniek van Hoek
%
% Goal:     Clusterings of US responses using k-means, performed on
%           co-registered neurons per session

% Input:    Change parameters below         
%           Saved clustering matrix can be loaded
%
% Output:   Heatmaps of co-registered neurons
%           Alluvial Plot
%           Traces per cluster

close all; 
clearvars -except clusters cluster_pre cluster_post cluster_reh

addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
cd = "Reesha/MATLAB/cohort_9_calcium_imaging/Recent";
%% set parameters
pre_event = 60; %time to take before bout start in seconds
post_event = 100; %set number to frames postlick to use to determine if the cell was responsive
z_score_window = 30; %number of frames in the extracted neural actiity window to use for zscore normalization
mice = ["M4-2", "M4-3", "M5-1", "M5-2", "M5-3"];
rerun = false; % Set to true if want to rerun clustering
sucrose = true;

%% Loading co registration of all mice
if sucrose
    dates = ["01082022", "02112022", "03012022"];
    reg_idx = cell(5,1);
    for mouse = 1:5
        filename = sprintf("celleg/Suc/%s_task_suc_cellreg.mat", mice(mouse));
        temp = load(fullfile(cd,filename));
        reg_idx{mouse} =  temp.cell_registered_struct.cell_to_index_map;
        reg_idx{mouse}(any(reg_idx{mouse}==0,2),:) = [];
        n_neurons_mouse(mouse) = size(reg_idx{mouse}(:,1),1);
    end
else
    dates = ["01232022", "02082022", "02282022"];
    reg_idx = cell(5,1);
    for mouse = 1:5
        filename = sprintf("celleg/Etoh/%s_task_etoh_cellreg.mat", mice(mouse));
        temp = load(fullfile(cd,filename));
        reg_idx{mouse} =  temp.cell_registered_struct.cell_to_index_map;
        reg_idx{mouse}(any(reg_idx{mouse}==0,2),:) = [];
        n_neurons_mouse(mouse) = size(reg_idx{mouse}(:,1),1);
    end
end

%% Loading neurons of all mice
trial_match = 10;
ses = cell(5,3);

for session = 1:3
    for mouse=1:5
      registration = reg_idx{mouse};

      name_zscores = sprintf('zscored/%s_%s_zscored_etoh.mat', dates(session), mice(mouse));
      temp = load(fullfile(cd, name_zscores));
      ses{mouse,session}   = temp.neuron_event_zscored;
      
      % Keep only cells coregistered over all 3  sessions
          for i =1:size(registration,1)
            neural_activity{mouse,session}(i,:,:) = ses{mouse,session}(registration(i,session),:,:);
            
            % Find mean trace over trial-matched trials
            for k = 1:size(neural_activity{mouse,session},3)
                mean_traces{mouse,session}(i,k) = mean(neural_activity{mouse,session}(i,1:trial_match,k));
              
            end
            
          end
      
    end
    full_traces{session} = vertcat(mean_traces{:,session});
end

%% Vector of length neurons with mouse ID
cell_ID = [];
for i = 1:size(neural_activity,1)
    cell_num = size(neural_activity{i,1}, 1);
    cell_ID = vertcat(cell_ID, i.*ones(cell_num,1));
end


%% Delete cell 39 and 5
if sucrose == false
    full_traces{1}(39,:)= [];
    full_traces{2}(39,:)= [];
    full_traces{3}(39,:)= [];

    full_traces{1}(5,:)= [];
    full_traces{2}(5,:)= [];
    full_traces{3}(5,:)= [];
    
    cell_ID(39,:) = [];
    cell_ID(5,:) = [];
end

%% K means clustering per session
if rerun
    cluster_pre = kmeans(full_traces{1}(:,[60:90]),4, 'Replicates', 20); % Set replicates to make sure it seems more deterministic (e.g. not other clusters everytime)
    cluster_post = kmeans(full_traces{2}(:,[60:90]),4,'Replicates', 20);
    cluster_reh = kmeans(full_traces{3}(:,[60:90]),4,'Replicates', 20);
    clusters = [cluster_pre, cluster_post, cluster_reh];

    [clusters_sort, idx] = sort(clusters);

    for session=1:3
    full_traces2{session} = full_traces{session}(idx(:,session),:);
    end

    figure
    Sessions = ["Pre isolation", "Post Isolation", "Regroup"];
    for session = 1:3
        subplot(1,3,session)
        imagesc(full_traces2{session}(:,1:90))
        xline(pre_event,'w','Linestyle','-','LineWidth',2)
        xline(z_score_window,'w','Linestyle','-','LineWidth',2)
        colormap(redbluecmap)
        title(Sessions(session));
        caxis([-5 5])
        colorbar
    end
% 
%% Alluvial plot
% Change order of clusters so numbers have the same meaning --> change this
% every time depending on the figure I get! Necessary for Alluvial plot to
% work.

    % If made a mistake, restore clusters variable before continuining
    clusters = [cluster_pre, cluster_post, cluster_reh];

    %Pre iso
    clusters_new(clusters(:,1)==1,1) = 4;
    clusters_new(clusters(:,1)==2,1) = 3;
    clusters_new(clusters(:,1)==3,1) = 2;
    clusters_new(clusters(:,1)==4,1) = 1;

    %clusters_new(clusters_new(:,1)==5,1) = 2;

    % Post iso
    clusters_new(clusters(:,2)==1,2) = 4;
    clusters_new(clusters(:,2)==2,2) = 3;
    clusters_new(clusters(:,2)==3,2) = 2;
    clusters_new(clusters(:,2)==4,2) = 1;

    %clusters_new(clusters_new(:,2)==5,2) = 2;

    % Rehousing
    clusters_new(clusters(:,3)==1,3) = 4;
    clusters_new(clusters(:,3)==2,3) = 3;
    clusters_new(clusters(:,3)==3,3) = 2;
    clusters_new(clusters(:,3)==4,3) = 1;
    %clusters_new(clusters(:,3)==4,3) = 2;
end


%% Continue script from here after loading data
if rerun == false
    if sucrose
        load("Aniek/Thesis Figures/Alluvial/clusters_suc_kmeans.mat")
    else
        load("Aniek/Thesis Figures/Alluvial/clusters_etoh_kmeans.mat")
    end

    [clusters_sort_new, idx_new] = sort(clusters_new);

    for session=1:3
    full_traces_new{session} = full_traces{session}(idx_new(:,session),:);
    end

    figure
    Sessions = ["Pre-isolation", "Isolation", "Post-isolation"];
    for session = 1:3
        subplot(1,3,session)
        imagesc(full_traces_new{session}(:,1:90))
        xline(pre_event,'w','Linestyle','-','LineWidth',2)
        xline(z_score_window,'w','Linestyle','-','LineWidth',2)
        colormap(redbluecmap)
        title(Sessions(session));
        caxis([-5 5])
        colorbar
    end

    % Create alluvial plot
    CreateSankeyPlot(clusters_new)
    
    % Saves which mouse clusters belong too   
    cluster_per_mouse = zeros(5,4,3);
    for session =1:3
        for mouse = 1:size(mice,2)
            for cluster = 1:4
                for i = 1:size(cell_ID,1)
                    if cell_ID(i,1) == mouse &&  clusters_new(i,session) == cluster
                        cluster_per_mouse(mouse,cluster,session) = cluster_per_mouse(mouse,cluster) + 1;
                    end
                end
            end
        end
    end
end

%% Mean and STD of cluster traces (for significance prism)

for session = 1:3
    for clusters = 1:4
        idx = find(clusters_sort_new(:,session)==clusters);
        for j = 1:length(idx)       
            data(j,:) = full_traces_new{session}(idx(j),:); 
        end
        data_mean{clusters}(:,session) = mean(data);
        data_sem{clusters}(:,session)  = std(data)/sqrt(size(data,1));
    end
end

        
        
%% Get traces out of clusters
if rerun == true
    full_traces_new = full_traces2;
    clusters_sort_new = clusters_sort;
end

label = [];
data = [];

% The size of the matrix
M = 3;
N = 4;
% The index according to your preferred ordering (column-wise)
i_colwise = 4;
% Conversion function
[jj,ii] = ind2sub([N,M],i_colwise);
i_rowwise = sub2ind([M,N],ii,jj); 

counter = 0;
figure
for i = 1:4
for session = 1:3
    
   % for i = 1:4 % number of clusters
        counter = counter + 1;
        idx = find(clusters_sort_new(:,session)==i);

        for j = 1:length(idx)       
            data(j,:) = full_traces_new{session}(idx(j),:);   
        end

        err = std(data)/sqrt(size(data,1));

        subplot(4,3, counter);
        [hl, he] = errorbar_pn(1:size(full_traces_new{session},2), mean(data), err, 123, 0.4);

        hold on
        cluster = num2str(i);
       % title(cluster);
        xline(30,'k','Linestyle','--','LineWidth',1);
        xline(60, 'r','Linestyle','--','LineWidth',1);
        yline(0, 'k', 'Linestyle', '--', 'LineWidth',1);
        xlim([0 91]);
      %  ylim([-2.5 5]);
      %  clear('data')
        xticks([]);
       % if i ~=1
          %  yticks({});
       % end
    end
end

%% Analyses below where to statistically compare transitions between sessions, not used in thesis
% Turnover matrix observed
ensemble_matrix_3 = zeros(4,4,4);
for i=1:size(clusters_new,1)
    for m = 1:4
        for n = 1:4
            for k = 1:4
                if clusters_new(i,1) == m && clusters_new(i,2) == n && clusters_new(i,3) == k
                    ensemble_matrix_3(m,n,k) = ensemble_matrix_3(m,n,k) + 1;
                end
            end
        end
    end
end

% Turnover matrix expected 
pre(1) = sum(clusters_new(:,1)==1);
pre(2)= sum(clusters_new(:,1)==2);
pre(3) = sum(clusters_new(:,1)==3);
pre(4) = sum(clusters_new(:,1)==4);

for i =1:4
    transition_matrix(i,:,:) = ensemble_matrix_3(i,:,:)/sum(clusters_new(:,1)==i);
end

% Chance matrix 
chance_levels(1) = sum(clusters_new(:,1)==1)/16;
chance_levels(2) = sum(clusters_new(:,1)==2)/16;
chance_levels(3) = sum(clusters_new(:,1)==3)/16;
chance_levels(4) = sum(clusters_new(:,1)==4)/16;

for i = 1:4
    chance_matrix(i,:,:) = repmat(chance_levels(i),4,4);
end

% Statistics 
aa3 = log(ensemble_matrix_3./chance_matrix);
aa3(isinf( aa3)) =0; % zero out infinite values
aa3(isnan( aa3)) =0; % zero out infinite values
df=63;
gstat2 = 2* sum(sum(sum(ensemble_matrix_3.*aa3))); % this is G^2 statistic
p2 = 1 - chi2cdf(gstat2,df); 




%% Building expected matrix using sequential log likelihood analyses
for i=1:4
    for j=1:4
        AB(i,j) = sum(x(i,j,:));
    end
end

for i=1:4
    for j=1:4
        BC(j,i) = sum(x(:,i,j));
    end
end

for i=1:4
    for j=1:4
        AC(j,i) = sum(x(i,:,j));
    end
end

E12 = [4 4 4 4; 4 4 4 4; 4 4 4 4; 4 4 4 4];
E13 = AB ./ E12;

E22 = [sum(E13,1);sum(E13,1);sum(E13,1);sum(E13,1)];
E23 = BC ./E22;

E23a = E13.*E23(1,:);
E23b = E13.*E23(2,:);
E23c = E13.*E23(3,:);
E23d = E13.*E23(4,:);

E32a = sum(E23a,2);
E32b = sum(E23b,2);
E32c = sum(E23c,2);
E32d = sum(E23d,2);

E32 = [E32a E32b E32c E32d];
E33 = AC.' ./E32;

Efinal(:,:,1) = E23a.*E33(:,1);
Efinal(:,:,2) = E23b.*E33(:,2);
Efinal(:,:,3) = E23c.*E33(:,3);
Efinal(:,:,4) = E23d.*E33(:,4);


%% Chi squared statistics on chain
x = ensemble_matrix_3;
aa3 = log(x./Efinal);
aa3(isinf( aa3)) =0; % zero out infinite values
aa3(isnan( aa3)) =0; % zero out infinite values
df=63;
gstat2 = 2* sum(x.*aa3); % this is G^2 statistic
p2 = 1 - chi2cdf(gstat2,df); 