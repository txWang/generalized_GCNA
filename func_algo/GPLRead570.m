% read probeID and gene symbol of a platform
function [probeID, geneSymbol, geneBankAcc] = GPLRead_570(annotationFile)

fid = fopen(annotationFile);

tmpL = [];
%%%%%%%%% modify here based on the files %%%%%%%%%%%%
while (length(strfind(tmpL, '#Gene Ontology Molecular Function')) == 0) &(feof(fid)~=1)
    tmpL = fgetl(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% annotation table begins
tmpL = fgetl(fid);
headLine = textscan(tmpL, '%s', 'delimiter', '\t');
headLine = headLine{1};
nColumns = length(headLine);
textStr = strcat(repmat('%s ', 1, nColumns));
A = textscan(fid, textStr, 'delimiter', '\t');
fclose(fid);

probeID = A{1};

%%%%%%%%% modify here based on the files %%%%%%%%%%%%
ind1 = strmatch('Gene Symbol', headLine);
ind2 = strmatch('GB_ACC', headLine);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GPL570 has multiple Gene Symbol corresponding to a probe
geneSymbol = A{ind1};
for i=1:length(geneSymbol)
    tmp = strsplit(geneSymbol{i},' /// ');
    geneSymbol{i} = tmp{1};
end
geneBankAcc = A{ind2};
end