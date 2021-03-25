function [twoSigma,conf95] = getSigmaConf(coef)

meanCoef = mean(coef);
stdCoef = std(coef);
n = length(coef);
sortedCoef = sort(coef);

[lowLimit, uppLimit] = deal(sortedCoef(floor(0.025*n)),sortedCoef(ceil(0.975*n)));  

twoSigma = [meanCoef - stdCoef, meanCoef + stdCoef];
conf95 = [lowLimit, uppLimit];
end

