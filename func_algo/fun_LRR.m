function cMatrix = fun_LRR(finalExp, lambda)

finalExp = normr(finalExp);
[coefMat,E2] = solve_lrr(finalExp',finalExp',lambda,1,1);
CAbs = abs(coefMat);
[Srt,Ind] = sort( CAbs,1,'descend' );
for i = 1:size(CAbs,2)
    CAbs(:,i) = CAbs(:,i) ./ (CAbs(Ind(1,i),i)+eps);
end
cMatrix = CAbs + CAbs';
end
