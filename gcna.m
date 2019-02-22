clc
clear

folder = 'data';
geoID = 'GSE31684';
str = load(fullfile(folder, geoID,[geoID,'.mat']));
finalExp = str.finalExp;
finalSym = str.finalSym;

% Size of clusters for futher analysis
minSize = 10; maxSize = 500;

%%%%%%%%%%%%%%%%%
% LRR based lmQCM
fprintf('LRR based lmQCM \n');
saveFolder = fullfile(folder, geoID, 'LRR'); mkdir(saveFolder);
% parameters
gamma = 0.2; lambda = 0.01; minClusterSize = 5;

cMatrix = fun_LRR(finalExp, lambda);
mergedCluster = fun_lmQCM(cMatrix, gamma, minClusterSize);
mergedCluster = mergedCluster(cellfun(@length, mergedCluster) >= minSize & cellfun(@length, mergedCluster) <= maxSize);
mergedCluster = cellfun(@(x) finalSym(x), mergedCluster, 'UniformOutput', false);
save(fullfile(saveFolder, 'cls.mat'),'mergedCluster');
fprintf('Number of co-expression modules identified: %d \n', length(mergedCluster));

%%%%%%%%%%%%%%%%%
% PCC based lmQCM
fprintf('PCC based lmQCM \n');
saveFolder = fullfile(folder, geoID, 'PCC'); mkdir(saveFolder);
% parameters
gamma = 0.2; minClusterSize = 5;

cMatrix = PCCMatrix(finalExp);
mergedCluster = fun_lmQCM(cMatrix, gamma, minClusterSize);
mergedCluster = mergedCluster(cellfun(@length, mergedCluster) >= minSize & cellfun(@length, mergedCluster) <= maxSize);
mergedCluster = cellfun(@(x) finalSym(x), mergedCluster, 'UniformOutput', false);
save(fullfile(saveFolder, 'cls.mat'),'mergedCluster');
fprintf('Number of co-expression modules identified: %d \n', length(mergedCluster));