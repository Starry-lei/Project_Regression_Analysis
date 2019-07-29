%-------------------------------------------------------------------------
% Find RSquared for a given model function and measured data
%
% Input: mdata - Measured data as a nxm-matrix
% fdata - Values of the model function f(t,x) with computed
% parameter x1,...,xn and the k columns of mdata
% Output: r2
%--------------------------------------------------------------------------
function r2 = RSquared(mdata,fdata)
[Q,P] = size(mdata);
yq = mean(mdata(:,P));
sf = sum((fdata - yq).^2);
mf = sum((mdata(:,P) - yq).^2);
r2 = sf/mf;