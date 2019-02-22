function cMatrix = PCCMatrix(finalExp)
%%%%%%%%%%%%%%%%%% Step 1 - Compute PCC matrix %%%%%%%%%%%%%%%%%%%%%%
if sum(sum(isnan(finalExp))) == 0 % no NaN in finalExp
    cMatrix = massivePCCWithoutNaN(finalExp);
else
    disp('NaN in finalExp\n');
    return;
end
cMatrix(1 : size(cMatrix, 1) + 1 : end) = 0; %remove self correlation

%%%% normalize PCC %%%%%
cMatrix = abs(cMatrix);
D = sum(cMatrix);
D_half = 1./sqrt (D);
for i = 1 : size(cMatrix, 1)
    cMatrix(i, :) = cMatrix(i, :) * D_half(i);
end
for i = 1 : size(cMatrix, 1)
    cMatrix(:, i) = cMatrix(:, i) * D_half(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
