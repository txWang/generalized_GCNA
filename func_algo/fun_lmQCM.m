function mergedCluster = fun_lmQCM(cMatrix, gamma, minClusterSize)

%parameters
%%% lmQCM
lambda = 1; t = 1;
%%% merge highly overlapped sub-networks
beta = 0.4; minClusterSizes = 5;

C = localMaximumQCM(abs(cMatrix), gamma, t, lambda);
for i = 1:length(C)
    C{i} = unique(C{i});
end
mergedCluster = mergeOverlapCls(C, beta, minClusterSize);

end