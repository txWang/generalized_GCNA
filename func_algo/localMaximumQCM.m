%%% Local maximal Quasi-Clique Merger %%%%
function C = localMaximumQCM(cMatrix, gamma, t, lambda)

C = {};
% find the local maximal edges
% cMatrix(i,j), correlation between gene i and gene j
% maxInd(i), strongest correlated gene to i
[maxV, maxInd] = max(cMatrix);
maxEdges = []; maxW = [];

% maxInd(i) = j
% j is strongest correlated gene to i
% maxV(i)== max(j, :)
% i is strongest correlated gene to j
for i = 1 : length(maxInd)
    if maxV(i)== max(cMatrix(maxInd(i), :))
        maxEdges = [maxEdges; maxInd(i), i];
        maxW = [maxW; maxV(i)];
    end
end

[sortMaxV, sortMaxInd] = sort(maxW, 'descend');
sortMaxEdges = maxEdges(sortMaxInd, :);

% fprintf('Local maximum edges: %d \n',length(sortMaxInd));

currentInit = 1; noNewInit = 0;
nodesInCluster = [];

% iteratively grow the cluster
while (currentInit <= length(sortMaxInd)) & (noNewInit == 0)
    % only consider correlation that is larger than gamma percent of 
    % the strongest correlation
    if (sortMaxV(currentInit) < gamma * sortMaxV(1))
        noNewInit = 1;
    else
        %both node on edge not in current clusters
        if (ismember(sortMaxEdges(currentInit, 1), nodesInCluster)==0 & ...
                ismember(sortMaxEdges(currentInit, 2), nodesInCluster)==0)
            newCluster = sortMaxEdges(currentInit, :); % two nodes
            addingMode = 1;
            currentDensity = sortMaxV(currentInit);
            nCp = 2;
            totalInd = 1 : size(cMatrix, 1);
            remainInd = setdiff(totalInd, newCluster);
            while addingMode == 1
                % sum of correlation between i and nodes in newCluster 
                neighborWeights = sum(cMatrix(newCluster, remainInd));
                
                [maxNeighborWeight, maxNeighborInd] = max(neighborWeights);
                % contribute(v,C)
                % the ratio of the edge weight increase of G(C) on adding the vertex v
                c_v = maxNeighborWeight/nCp; 
                alphaN = 1 - 1/(2*lambda*(nCp+t));
                % threshold = alphaN * density (G(C))
                if (c_v >= alphaN * currentDensity)
                    % add v to newCluster if contribute(v,C) >= threshold
                    newCluster = [newCluster, remainInd(maxNeighborInd)];
                    nCp = nCp+1;
                    % fully connected undirected graph, v*(v-1)/2 edges
                    currentDensity = (currentDensity*((nCp-1)*(nCp-2)/2)+maxNeighborWeight)/(nCp*(nCp-1)/2);
                    remainInd = setdiff(remainInd, remainInd(maxNeighborInd));
                else
                    addingMode = 0;
                end
            end
            nodesInCluster = [nodesInCluster, newCluster];
            C = [C, newCluster];
            
        end
    end
    currentInit = currentInit + 1;
end
% fprintf('Cluster center considered: %d \n', currentInit-1);
end