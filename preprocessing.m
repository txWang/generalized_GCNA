clc;
clear;

topN = 10000; %number of unique genes to retain

folder = 'data';
geoID = 'GSE31684';
geoFile = 'GSE31684_series_matrix.txt';
GPLFile = 'GPL570-55999.txt';

[probeID, geneSymbol, geneBankAcc] = GPLRead570(fullfile(folder, geoID, GPLFile));
geoData = geoseriesread(fullfile(folder, geoID,geoFile));

%%%% match probe with gene symbol %%%
[probeID3,idx1,idx2] = intersect(probeID,geoData.Data.rownames,'legacy');
geneSym = geneSymbol(idx1);
geneData = geoData.Data(idx2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[finalExp, finalSym, uniExp, uniSym] = dataPreProcessing(geneData, geneSym, topN);
save(fullfile(folder, geoID,[geoID,'.mat']),'finalExp','finalSym','uniExp','uniSym','-v7.3');
