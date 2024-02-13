function dxdt = differential_equation(t,x,param)

% this differential equation is not an explicit function of time, but
% Matlab will later require that to be there for other integration
% functions

% This function returns the rate of change of x as a function of x and t
% and the parameters of the equation

dxdt = param.r ...
    + param.beta*x.^param.n./(param.x_M.^param.n + x.^param.n) ...
    - param.alpha*x;