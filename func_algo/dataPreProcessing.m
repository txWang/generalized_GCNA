%%%%%%%%%%%%%%%%%% Pre Processing %%%%%%%%%%%%%%%%%%%%%%
% filter out low mean and low variance gene based on ALL data
% default
function [finalExp, finalSym, uniExp, uniSym] ...
    = dataPreProcessing(geoDataSort, geneSymSort, topN)

if mean(mean(geoDataSort)) > 30 %need to take log
    fprintf('Take log2 of the orginal data.\n');
    geoDataSort = log2(geoDataSort);
end
%remove multiple gene exp values, retain only the gene exp value with
%highest mean for each gene using code HighExpressionProbes.m
[indUni, uniSym] = HighExpressionProbes(geneSymSort, geneSymSort, geoDataSort);
uniExp = double(geoDataSort);
uniExp = uniExp(indUni,:);
%delete empty symbol
idxEmpty = find(cellfun(@isempty,uniSym));
uniSym(idxEmpty) = [];
uniExp(idxEmpty,:) = [];
%This is background

%remove data with lowest 20% absolute exp value shared by all samples
[mask, geoDataFilter, geneSymFilter]= genelowvalfilter(geoDataSort,geneSymSort,'percentile',20);
%remove data with lowest 10% variance across samples
[mask, geoDataFilter2, geneSymFilter2] = genevarfilter(geoDataFilter,geneSymFilter);
expData = double(geoDataFilter2);
%remove multiple gene exp values, retain only the gene exp value with
%highest mean for each gene using code HighExpressionProbes.m
[ind1, tmpSym] = HighExpressionProbes(geneSymFilter2, geneSymFilter2, expData);
tmpExp = expData(ind1,:); %rows re-sorted to the alphabetic order of uniGene.

%sort by mean expression
[sortMean, sortInd] = sort(mean(tmpExp, 2), 'descend');
finalExp = tmpExp(sortInd, :); 
finalSym = tmpSym(sortInd);

%delete empty symbol
idxEmpty = find(cellfun(@isempty,finalSym));
finalSym(idxEmpty) = [];
finalExp(idxEmpty,:) = [];

%choose genes with higest topN expression
finalExp = finalExp(1:topN, :); 
finalSym = finalSym(1:topN);

fprintf('NaN in final expression: %d\nNaN in background expression: %d\n',...
    sum(sum(isnan(finalExp))), sum(sum(isnan(uniExp))));
end