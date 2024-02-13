function dxdt = differential_equation_timedep(t,x,param)

% this differential equation is not an explicit function of time, but
% Matlab will later require that to be there for other integration
% functions

% This function returns the rate fo change of x as a function of x and t
% and the parameters of the equation

% Please complete this line of code to step rt (the local param.r_basal)
% using the function 'steppingfunction'
rt = param.r*steppingfunction(t); % remove the '0' and put in a nice function here.

dxdt = rt ...
    + param.beta*x.^param.n./(param.x_M.^param.n + x.^param.n) ...
    - param.alpha*x;