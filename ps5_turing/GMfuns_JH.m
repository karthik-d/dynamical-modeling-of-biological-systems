<<<<<<< HEAD
%%    GMfuns_JH
=======
%    GMfuns_JH
>>>>>>> c988670df76dcf42daf1b702dd605bd98907b7a9
%   
%in this script, we define the necessary functions for our Gierer-Meinhardt
%system:
%1) the initial conditions;
%2) the boundary conditions;
%3) the PDEs themseleves.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here we use function handles to define subfunctions. This makes life
%easier. 
function funs = GMfuns_JH
    funs.fun1 = @GMICfun;
    funs.fun2 = @GMBCfun;
    funs.fun3 = @GMPDEfun;
end

%initial conditions function
function u0 = GMICfun(x,P)
% Here we define the initial conditions for the system --
% We use a small random perturbation from the steady state for
%the activator and do not perturb the inhibitor (note: steady state
%concentrations are 1 for both activator & inhibitor).
u0 = [ones(1,length(x)) +  rand(1,length(x))/100 ; ones(1,length(x)) ];

end

%boundary conditions function
function [pl,ql,pr,qr] = GMBCfun(xl,ul,xr,ur,t,P)
% Here we define the boundary conditions for the system (no flux at either end)--
pl = [0;0]; ql = [1;1];
pr = [0;0]; qr = [1;1];
end

%PDE
function [c,f,s] = GMPDEfun(x,t,u,dudx,P)
%Here we define the PDE
D = P(1);
omega1 = P(2);
sigm = P(3);
% PDE
c = [1;1];
f = [D;1].*dudx;
s = [ u(1)^2/(0.2+u(2)) - u(1) + sigm ; omega1*(u(1)^2 - u(2)) ];
end