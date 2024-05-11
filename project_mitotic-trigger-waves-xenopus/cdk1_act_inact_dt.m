function [pde_funcs] = cdk1_act_inact_dt()
    pde_funcs.Ckd1_PDE = @Cdk1_PDE;
    pde_funcs.Ckd1_init = @Cdk1_init;
    pde_funcs.Ckd1_boundaries = @Cdk1_boundaries;
end


% PDE Function.
function [c,f,s] = Cdk1_PDE(x, t, u, dudx, consts)

    % ** Independent vars are encoded as: u = [Cdk1_act; Cdk1_inact] **
        
    % === Define dependent variables. ===
    % diffusion constant in micro-m^2 min^-1.
    consts.Dp = 600;

    % syntesis const.
    ksynth = 1.5;

    % rate constants in min^-1.
    % -  for Cdc25.
    % draw from x-dependent function.
    acdc25 = consts.scale * acdc25_at_x(x);
    bcdc25 = consts.scale * bcdc25_at_x(x);
    
    % -  for Wee1.
    awee1 = consts.scale * 0.4;
    bwee1 = consts.scale * 2;
    
    % -  for degradation.
    adeg = 0.01;
    bdeg = 0.06;
    
    % Hill coefficients.
    ncdc25 = 11;
    nwee1 = 3.5;
    ndeg = 17;

    % === Define PDE coefficients. ===
    
    % time coupling term.
    c = [1; 1]

    % flux term.
    f = [consts.Dp*dudx(1); consts.Dp*dudx(2)];

    % source term.
    s = [   ... 
			ksynth ...
			+ (acdc25 + bcdc25*(u(1)^ncdc25)/(consts.EC50cdc25^ncdc25 + u(1)^ncdc25))*u(2) ...
			- (awee1 + bwee1*(consts.EC50wee1^nwee1)/(consts.EC50wee1^nwee1 + u(1)^nwee1))*u(1) ...
			- (adeg + bdeg*(u(1)^ndeg)/(consts.EC50deg^ndeg + u(1)^ndeg))*u(1) ...
		;
			- (acdc25 + bcdc25*(u(1)^ncdc25)/(consts.EC50cdc25^ncdc25 + u(1)^ncdc25))*u(2) ...
			+ (awee1 + bwee1*(consts.EC50wee1^nwee1)/(consts.EC50wee1^nwee1 + u(1)^nwee1))*u(1) ...
			- (adeg + bdeg*(u(1)^ndeg)/(consts.EC50deg^ndeg + u(1)^ndeg))*u(2) ...
        ];
end


function [u0] = Ckd1_init(x, consts)
    u0 = [0; 0];
end


function [pl, ql, pr, qr] = Ckd1_boundaries(xl, ul, xr, ur, t, consts)
    % boundary condition: flux is 0 at both ends - xl and xr.
    pl = 0;
    ql = 1/consts.Dp;
    % -- 
    pr = 0;
    qr = 1/consts.Dp;
end


% functions to encode the inhomogeneity induced by centrome
% as 50% higher concentration in the middle 5-micron section.
function [acdc25] = acdc25_at_x(x)
    if and(x>=-2.5, x<2.5)
        acdc25 = 0.12;
    else
        acdc25 = 0.08;
    end
end

function [bcdc25] = bcdc25_at_x(x)
    if and(x>=-2.5, x<2.5)
        bcdc25 = 0.6;
    else
        bcdc25 = 0.4;
    end
end


%{

=== Mathematica Implementation for Reference. ====

NDSolve[
{
    D[Cdk1[t, x], t] ==
        +ksynth[t, x]
        + (acdc25[t, x]
        + bcdc25[t, x] Cdk1[t, x]^11/(
        EC50cdc25^11 + Cdk1[t, x]^11)) (Cdk1inact[t, x])
        - (awee1[t, x]
        + bwee1[t, x] EC50wee1^3.5/(EC50wee1^3.5 + Cdk1[t, x]^3.5)) Cdk1[t, x]
        - (adeg[t, x]
        + bdeg[t, x] Cdk1[t, x]^17/(EC50deg^17 + Cdk1[t, x]^17)) Cdk1[t, x]
        + Dp D[Cdk1[t, x], x, x],
    
    D[Cdk1inact[t, x], t] ==
        - (adeg[t, x] +
        bdeg[t, x] Cdk1[t, x]^17/(EC50deg^17 + Cdk1[t, x]^17)) (Cdk1inact[t, x])
        - (acdc25[t, x] +
        bcdc25[t, x] Cdk1[t, x]^11/(
        EC50cdc25^11 + Cdk1[t, x]^11)) (Cdk1inact[t, x])
        + (awee1[t, x] +
        bwee1[t, x] EC50wee1^3.5/(EC50wee1^3.5 + Cdk1[t, x]^3.5)) Cdk1[t, x]
        + Dp D[Cdk1inact[t, x], x, x],
    
    D[adeg[t, x], t] == 0,
    D[bdeg[t, x], t] == 0,
    D[acdc25[t, x], t] == 0,
    D[bcdc25[t, x], t] == 0,
    D[awee1[t, x], t] == 0,
    D[bwee1[t, x], t] == 0,
    D[ksynth[t, x], t] == 0,
    
    ksynth[0, x] == 1.5(*2.5*),
    adeg[0, x] == .01,
    bdeg[0, x] == .06,
    awee1[0, x] == scale*.08*.5,
    bwee1[0, x] == scale*.4*.5,
    acdc25[0, x] == scale*Piecewise[
    {
        {0.08, x < -2.5},
        {.12, -2.5 <= x < 2.5},
        {0.08, 2.5 <= x}
    }],
    bcdc25[0, x] == scale*Piecewise[
    {
        {0.4, x < -2.5},
        {.6, -2.5 <= x < 2.5},
        {0.4, 2.5 <= x}
    }],
    Cdk1inact[0, x] == 0,
    Cdk1[0, x] == 0,
    (D[Cdk1[t, x], x] /. x -> -xrange) == 0,
    (D[Cdk1[t, x], x] /. x -> xrange) == 0,
    (D[Cdk1inact[t, x], x] /. x -> -xrange) == 0,
    (D[Cdk1inact[t, x], x] /. x -> xrange) == 0
}, --> EQUATIONS SYSTEM.
{
    Cdk1, Cdk1inact, acdc25, bcdc25, awee1, bwee1, ksynth, adeg, bdeg
}, --> DEPENDENT VARS. 
{
t, 0, tmax
}, --> TIME INTERVAL.
{
    x, -xrange, xrange
}]; --> INDEPENDENT VARS.

%}