%% Plots of pseudopopulation decodings
% Aniek van Hoek
%
% Goal:     Script gets information out of saved decoding cell arrays and
%           plots AUC and accuracy per model, comparing test and shuffled data

% Input:    Manually load the decoding results in
%           Reesha/MATLAB/cohort_9_calcium_imaging/Recent/Decoding
%           Change data variable depending on which model performance to retrieve
%
% Output:   Boxplots of decoding performance for either generalization over
%           sessions or for within session decoding
%           sig = significance using ranksum test

clearvars; clc; close all;
addpath(genpath(fullfile('/nadata', 'snlkt', 'data', 'Aniek', 'Code')));
%% Creating new cell arrays that contain the decoding performance
%% Ethanol
counter = 0;
for data = [1,3]
    counter = counter+1;
    etoh_AUC_ses1{1,1}(:,counter) = result_watervsetoh{data, 1};
    etoh_AUC_ses1{1,2}(:,counter) = result_watervsetoh_shuffle{data, 1};
    sig_AUC_etoh_ses1(counter) = ranksum(etoh_AUC_ses1{1,1}(:,counter), etoh_AUC_ses1{1,2}(:,counter));

    etoh_AUC_ses2{1,1}(:,counter) = result_watervsetoh{data, 2};
    etoh_AUC_ses2{1,2}(:,counter) = result_watervsetoh_shuffle{data, 2};
    sig_AUC_etoh_ses2(counter) = ranksum(etoh_AUC_ses2{1,1}(:,counter), etoh_AUC_ses2{1,2}(:,counter));
    
    etoh_AUC_ses3{1,1}(:,counter) = result_watervsetoh{data, 3};
    etoh_AUC_ses3{1,2}(:,counter) = result_watervsetoh_shuffle{data, 3};
    sig_AUC_etoh_ses3(counter) = ranksum(etoh_AUC_ses3{1,1}(:,counter), etoh_AUC_ses3{1,2}(:,counter));
    
    etoh_sig_1vs2_auc(counter) = ranksum(result_watervsetoh{data, 1},  result_watervsetoh{data, 2});
    etoh_sig_2vs3_auc(counter) = ranksum(result_watervsetoh{data, 2},  result_watervsetoh{data, 3});
    etoh_sig_1vs3_auc(counter) = ranksum(result_watervsetoh{data, 1},  result_watervsetoh{data, 3});
end

counter = 0;
for data = [2,4]
    counter = counter+1;
    etoh_ACC_ses1{1,1}(:,counter) = result_watervsetoh{data, 1}; 
    etoh_ACC_ses1{1,2}(:,counter) = result_watervsetoh_shuffle{data, 1};
    sig_ACC_etoh_ses1(counter) = ranksum(etoh_ACC_ses1{1,1}(:,counter), etoh_ACC_ses1{1,2}(:,counter));
    
    etoh_ACC_ses2{1,1}(:,counter) = result_watervsetoh{data, 2}; 
    etoh_ACC_ses2{1,2}(:,counter) = result_watervsetoh_shuffle{data, 2};
    sig_ACC_etoh_ses2(counter) = ranksum(etoh_ACC_ses2{1,1}(:,counter), etoh_ACC_ses2{1,2}(:,counter));

    etoh_ACC_ses3{1,1}(:,counter) = result_watervsetoh{data, 3}; 
    etoh_ACC_ses3{1,2}(:,counter) = result_watervsetoh_shuffle{data, 3};
    sig_ACC_etoh_ses3(counter) = ranksum(etoh_ACC_ses3{1,1}(:,counter), etoh_ACC_ses3{1,2}(:,counter));
    
    etoh_sig_1vs2_acc(counter) = ranksum(result_watervsetoh{data, 1},  result_watervsetoh{data, 2});
    etoh_sig_2vs3_acc(counter) = ranksum(result_watervsetoh{data, 2},  result_watervsetoh{data, 3});
    etoh_sig_1vs3_acc(counter) = ranksum(result_watervsetoh{data, 1},  result_watervsetoh{data, 3});
end


%% Sucrose
counter = 0;
for data = [3,5,7]
    counter = counter+1;
    suc_AUC_ses1{1,1}(:,counter) = result_watervssucrose{data, 1};
    suc_AUC_ses1{1,2}(:,counter) = result_watervssucrose_shuffle{data, 1};
    sig_AUC_suc_ses1(counter) = ranksum(suc_AUC_ses1{1,1}(:,counter), suc_AUC_ses1{1,2}(:,counter));
    
    suc_AUC_ses2{1,1}(:,counter) = result_watervssucrose{data, 2};
    suc_AUC_ses2{1,2}(:,counter) = result_watervssucrose_shuffle{data, 2};
    sig_AUC_suc_ses2(counter) = ranksum(suc_AUC_ses2{1,1}(:,counter), suc_AUC_ses2{1,2}(:,counter));
    
    suc_AUC_ses3{1,1}(:,counter) = result_watervssucrose{data, 3};
    suc_AUC_ses3{1,2}(:,counter) = result_watervssucrose_shuffle{data, 3};
    sig_AUC_suc_ses3(counter) = ranksum(suc_AUC_ses3{1,1}(:,counter), suc_AUC_ses3{1,2}(:,counter));
    
    suc_sig_1vs2_auc(counter) = ranksum(result_watervssucrose{data, 1},  result_watervssucrose{data, 2});
    suc_sig_2vs3_auc(counter) = ranksum(result_watervssucrose{data, 2},  result_watervssucrose{data, 3});
    suc_sig_1vs3_auc(counter) = ranksum(result_watervssucrose{data, 1},  result_watervssucrose{data, 3});
end

counter = 0;
for data = [4,6,8]
    counter = counter+1;
    suc_ACC_ses1{1,1}(:,counter) = result_watervssucrose{data, 1}; 
    suc_ACC_ses1{1,2}(:,counter) = result_watervssucrose_shuffle{data, 1};
    sig_ACC_suc_ses1(counter) = ranksum(suc_ACC_ses1{1,1}(:,counter), suc_ACC_ses1{1,2}(:,counter));
    
    suc_ACC_ses2{1,1}(:,counter) = result_watervssucrose{data, 2}; 
    suc_ACC_ses2{1,2}(:,counter) = result_watervssucrose_shuffle{data, 2};
    sig_ACC_suc_ses2(counter) = ranksum(suc_ACC_ses2{1,1}(:,counter), suc_ACC_ses2{1,2}(:,counter));
    
    suc_ACC_ses3{1,1}(:,counter) = result_watervssucrose{data, 3}; 
    suc_ACC_ses3{1,2}(:,counter) = result_watervssucrose_shuffle{data, 3};
    sig_ACC_suc_ses3(counter) = ranksum(suc_ACC_ses3{1,1}(:,counter), suc_ACC_ses3{1,2}(:,counter));
    
    suc_sig_1vs2_acc(counter) = ranksum(result_watervssucrose{data, 1},  result_watervssucrose{data, 2});
    suc_sig_2vs3_acc(counter) = ranksum(result_watervssucrose{data, 2},  result_watervssucrose{data, 3});
    suc_sig_1vs3_acc(counter) = ranksum(result_watervssucrose{data, 1},  result_watervssucrose{data, 3});
end


%% Figures Ethanol
% AUC
c = {'SVM Linear', 'Random Forest'}; 
figure
subplot(2,3,1)
boxplotGroup(etoh_AUC_ses1, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle",}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("ethanol pre-isolation")     

subplot(2,3,2)
boxplotGroup(etoh_AUC_ses2, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("ethanol isolation")    

subplot(2,3,3)
boxplotGroup(etoh_AUC_ses3, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("ethanol post-isolation")    

% Accuracy
subplot(2,3,4)
boxplotGroup(etoh_ACC_ses1, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("ethanol pre-isolation")     

subplot(2,3,5)
boxplotGroup(etoh_ACC_ses2, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("ethanol isolation")    

subplot(2,3,6)
boxplotGroup(etoh_ACC_ses3, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("ethanol post-isolation")  

%% Figures Sucrose
% AUC
c = {'SVM Linear', 'Random Forest', 'Linear classifier'}; 
figure
subplot(2,3,1)
boxplotGroup(suc_AUC_ses1, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("sucrose pre-isolation")     

subplot(2,3,2)
boxplotGroup(suc_AUC_ses2, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("sucrose isolation")    

subplot(2,3,3)
boxplotGroup(suc_AUC_ses3, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("sucrose post-isolation")    


% Accuracy
subplot(2,3,4)
boxplotGroup(suc_ACC_ses1, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("sucrose pre-isolation")     

subplot(2,3,5)
boxplotGroup(suc_ACC_ses2, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("sucrose isolation")    

subplot(2,3,6)
boxplotGroup(suc_ACC_ses3, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("sucrose post-isolation")  
        
%% Per session decoding 
counter = 0;
for data = [1,3,5]
    counter = counter+1;
    AUCses1vs2{1,1}(:,counter) = result_ses1vs2{1,data};
    AUCses1vs2{1,2}(:,counter) = result_ses1vs2_shuffle{1,data};
    sig_AUC_etoh_ses1vs2(counter) = ranksum(AUCses1vs2{1,1}(:,counter), AUCses1vs2{1,2}(:,counter));

    AUCses1vs3{1,1}(:,counter) = result_ses1vs3{1,data};
    AUCses1vs3{1,2}(:,counter) = result_ses1vs3_shuffle{1,data};
    sig_AUC_etoh_ses1vs3(counter) = ranksum(AUCses1vs3{1,1}(:,counter), AUCses1vs3{1,2}(:,counter));
    
    AUCses2vs3{1,1}(:,counter) = result_ses2vs3{1,data};
    AUCses2vs3{1,2}(:,counter) = result_ses2vs3_shuffle{1,data};
    sig_AUC_etoh_ses2vs3(counter) = ranksum(AUCses2vs3{1,1}(:,counter), AUCses2vs3{1,2}(:,counter));
    
    AUCses2vs1{1,1}(:,counter) = result_ses2vs1{1,data};
    AUCses2vs1{1,2}(:,counter) = result_ses2vs1_shuffle{1,data};
    sig_AUC_etoh_ses2vs1(counter) = ranksum(AUCses2vs1{1,1}(:,counter), AUCses2vs1{1,2}(:,counter));
    
    AUCses3vs1{1,1}(:,counter) = result_ses3vs1{1,data};
    AUCses3vs1{1,2}(:,counter) = result_ses3vs1_shuffle{1,data};
    sig_AUC_etoh_ses3vs1(counter) = ranksum(AUCses3vs1{1,1}(:,counter), AUCses3vs1{1,2}(:,counter));
    
    AUCses3vs2{1,1}(:,counter) = result_ses3vs2{1,data};
    AUCses3vs2{1,2}(:,counter) = result_ses3vs2_shuffle{1,data};
    sig_AUC_etoh_ses3vs2(counter) = ranksum(AUCses3vs2{1,1}(:,counter), AUCses3vs2{1,2}(:,counter));
   
end

counter = 0;
for data = [2,4,6]
    counter = counter+1;
    ACCses1vs2{1,1}(:,counter) = result_ses1vs2{1,data};
    ACCses1vs2{1,2}(:,counter) = result_ses1vs2_shuffle{1,data};
    sig_ACC_etoh_ses1vs2(counter) = ranksum(ACCses1vs2{1,1}(:,counter), ACCses1vs2{1,2}(:,counter));

    ACCses1vs3{1,1}(:,counter) = result_ses1vs3{1,data};
    ACCses1vs3{1,2}(:,counter) = result_ses1vs3_shuffle{1,data};
    sig_ACC_etoh_ses1vs3(counter) = ranksum(ACCses1vs3{1,1}(:,counter), ACCses1vs3{1,2}(:,counter));
    
    ACCses2vs3{1,1}(:,counter) = result_ses2vs3{1,data};
    ACCses2vs3{1,2}(:,counter) = result_ses2vs3_shuffle{1,data};
    sig_ACC_etoh_ses2vs3(counter) = ranksum(ACCses2vs3{1,1}(:,counter), ACCses2vs3{1,2}(:,counter));
    
    ACCses2vs1{1,1}(:,counter) = result_ses2vs1{1,data};
    ACCses2vs1{1,2}(:,counter) = result_ses2vs1_shuffle{1,data};
    sig_ACC_etoh_ses2vs1(counter) = ranksum(ACCses2vs1{1,1}(:,counter), ACCses2vs1{1,2}(:,counter));
    
    ACCses3vs1{1,1}(:,counter) = result_ses3vs1{1,data};
    ACCses3vs1{1,2}(:,counter) = result_ses3vs1_shuffle{1,data};
    sig_ACC_etoh_ses3vs1(counter) = ranksum(ACCses3vs1{1,1}(:,counter), ACCses3vs1{1,2}(:,counter));
    
    
    ACCses3vs2{1,1}(:,counter) = result_ses3vs2{1,data};
    ACCses3vs2{1,2}(:,counter) = result_ses3vs2_shuffle{1,data};
    sig_ACC_etoh_ses3vs2(counter) = ranksum(ACCses3vs2{1,1}(:,counter), ACCses3vs2{1,2}(:,counter));
   
end

sig_ACC_etoh_1vs2VS2vs1 = ranksum(result_ses1vs2{1,2}, result_ses2vs1{1,2});
%% Plots
% AUC
c = {'SVM Linear', 'Random Forest', 'Linear classifier'}; 
figure
subplot(2,3,1)
boxplotGroup(AUCses1vs2, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("ethanol session 1 vs 2")     

subplot(2,3,4)
boxplotGroup(AUCses1vs3, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("ethanol session 1 vs 3")   

subplot(2,3,2)
boxplotGroup(AUCses2vs3, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("ethanol session 2 vs 3")  

subplot(2,3,5)
boxplotGroup(AUCses2vs1, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("ethanol session 2 vs 1")  

subplot(2,3,3)
boxplotGroup(AUCses3vs1, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("ethanol session 3 vs 1")  

subplot(2,3,6)
boxplotGroup(AUCses3vs2, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("AUC")
title("ethanol session 3 vs 2")  


% Accuracy

figure
subplot(2,3,1)
boxplotGroup(ACCses1vs2, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("ethanol session 1 vs 2")     

subplot(2,3,4)
boxplotGroup(ACCses1vs3, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("ethanol session 1 vs 3")   

subplot(2,3,2)
boxplotGroup(ACCses2vs3, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("ethanol session 2 vs 3")  

subplot(2,3,5)
boxplotGroup(ACCses2vs1, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("ethanol session 2 vs 1")  

subplot(2,3,3)
boxplotGroup(ACCses3vs1, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("ethanol session 3 vs 1")  

subplot(2,3,6)
boxplotGroup(ACCses3vs2, 'groupLines', true, 'interGroupSpace',1, 'primaryLabels', {"test", "shuffle","test", "shuffle","test", "shuffle"}, ...
    'secondaryLabels', c, ...
    'groupLabelType', 'horizontal')  
ylabel("ACC")
title("ethanol session 3 vs 2")  
