% Elif Ayvali 11/03/2015 eayvali@gmail.com
% Upper Confidence Bound
% ymu: the mean of GP posterior
% ys: the standard deviation of GP posterior
function res = UCB(ymu,ys2,k)

switch nargin
    case 2
        beta=1.96;
    case 3
        beta = k;
end

ys=sqrt(ys2);
res =ymu+beta.*ys;
end
