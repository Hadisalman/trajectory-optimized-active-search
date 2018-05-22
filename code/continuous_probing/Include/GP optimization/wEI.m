% Elif Ayvali 06/16/2015 eayvali@gmail.com
% yEI: the value of y to 'improve over'.
% ymu: the mean of GP posterior
% ys: the standard deviation of GP posterior
% w = 0 global exploration 
% w = 1 local exploitation
%w=0.5 wEI becomes EI
function res = wEI(yEI,ymu,ys2,w)
eps=0.01;
ys=sqrt(ys2);
res = w*(ymu-yEI-eps).*normcdf((ymu-yEI-eps)./ys)+(1-w)*ys.*normpdf((ymu-yEI-eps)./ys);
end
