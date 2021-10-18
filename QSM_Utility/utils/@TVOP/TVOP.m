function  res = TVOP()

%res = TVOP()
%
% Implements a spatial finite-differencing operator.
%

res.adjoint = 0;
res = class(res,'TVOP');

