%%%%%%%%%%%%%%%%%%%%%%
%Output the probe with highest mean among all samples of a gene
function [ind, uniGene] = HighExpressionProbes(Genes, Probes, Data)
meanData = mean(Data, 2);
[sGenes, sInd] = sort(Genes);
sProbes = Probes(sInd);
sMean = meanData(sInd);
% To use Matlab 2012b and before, the 'legacy' option allows the output to
% be the last occurence of each gene
% [uniGene, uniInd] = unique(sGenes, 'legacy');
% uniInd = [0; uniInd];

%%% Using the new version of the "unique()" function, which output the
%%% indices of the first occurence of the genes
[uniGene, uniInd] = unique(sGenes);
uniInd = [uniInd; length(sGenes)+1];

tmpInd = zeros(1, length(uniGene));
for i = 1 : length(uniInd)-1
    % for Matlab 2012b and before
    % [maxV, maxInd] = max(sMean(uniInd(i)+1:uniInd(i+1)));
    
    % for new versions of Matlab
    [maxV, maxInd] = max(sMean(uniInd(i):uniInd(i+1)-1));
    tmpInd(i) = uniInd(i) + maxInd - 1;
end
% [maxV, maxInd] = max(sMean(uniInd(end):end));
% tmpInd(end) = uniInd(end) + maxInd - 1;
ind = sInd(tmpInd);
end