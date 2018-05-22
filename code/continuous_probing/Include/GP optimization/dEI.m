% Elif Ayvali 06/16/2015 eayvali@gmail.com
% yEI: the value of y to 'improve over'.
% ymu: the mean of GP posterior
% ys: the standard deviation of GP posterior
function res = EI(yEI,ymu,ys2)
eps=0.01;
ys=sqrt(ys2);
res = (ymu-yEI-eps).*normcdf((ymu-yEI-eps)./ys)+ys.*normpdf((ymu-yEI-eps)./ys);
end
