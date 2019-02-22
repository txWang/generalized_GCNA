function mergedCluster = mergeOverlapCls(C, beta, minClusterSize)

%%%%%%%% Step 3 - Merge the overlapped networks %%%%%%%%%%%%%
% Merging highly overlapped sub-networks in C with respect to beta

%%%%%%%% Sort and iteratively merge %%%%%%%%%%%%%%%%%%%%%%%%%
sizeC = zeros(1, length(C));
for i = 1 : length(C)
    sizeC(i) = length(C{i});
end

[sortC, sortInd] = sort(sizeC, 'descend');
C = C(sortInd);

ind = find(sortC >= minClusterSize);

mergedCluster = C(ind);
mergeOccur = 1; 
currentInd = 0;

% if merged happend in this loop, start over again
while mergeOccur == 1
    mergeOccur = 0;
    while currentInd < length(mergedCluster)
        currentInd = currentInd + 1;
        excludeInd = [];
        if (currentInd < length(mergedCluster))
            % cluster that could be merged into 1 : currentInd has been
            % merged, consider clusters starting from currentInd+1
            keepInd = 1 : currentInd;
            % whether cluster j could be merged with cluster currentInd
            for j = currentInd+1 : length(mergedCluster)
                interCluster = intersect(mergedCluster{currentInd}, mergedCluster{j});
                % if proportion of intersect elements larger than threshold
                % beta, merge
                if length(interCluster) >= beta*min(length(mergedCluster{j}), length(mergedCluster{currentInd}))
                    mergedCluster{currentInd} = union(mergedCluster{currentInd}, mergedCluster{j});
                    mergeOccur = 1;
                else
                    keepInd = [keepInd, j];
                end;
            end;
            mergedCluster = mergedCluster(keepInd);
            %fprintf('current size of merged clusters: %d \n',length(mergedCluster));
        end;
    end;
    % sort mergedCluster according to size
    sizeMergedCluster = zeros(1, length(mergedCluster));
    for i = 1 : length(mergedCluster)
        sizeMergedCluster(i) = length(mergedCluster{i});
    end;
    [sortSize, sortMergedInd] = sort(sizeMergedCluster, 'descend');
    mergedCluster = mergedCluster(sortMergedInd);
    currentInd = 0;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end