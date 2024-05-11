function funcs = cdk1_cyclin_ode
    funcs.dydt = @cdk1_cyclin_ode_dydt;
    % funcs.nullclines = @cdk1_cyclin_ode_nullclines;
end

%{

=== Parameter Values ===
    ksynth	1 nM/min
    adeg	0.01 min-1
    bdeg	0.04 min-1
    EC50deg	32 nM
    ndeg	17
    aCdc25	0.16 min-1
    bCdc25	0.80 min-1
    EC50Cdc25	35 nM
    nCdc25	11
    aWee1	0.08 min-1
    bWee1	0.40 min-1
    EC50Wee1	30 nM
    nWee1	3.5

%}


function dydt = cdk1_cyclin_ode_dydt(t, y, params, consts)

    syms ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc;

    % dydt is a vector of rates for [Cdk1; Cyclin].
    [cdk1_expr, cyclin_expr] = cdk1_cyclin_ode_exprs;
    dydt = [ double(subs(cdk1_expr, {ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc}, ...
                {params.ksynth params.acdc25 params.bcdc25 params.ncdc25 params.awee1 params.bwee1 params.nwee1 params.adeg params.bdeg params.ndeg ...
                consts.EC50cdc25 consts.EC50wee1 consts.EC50deg y(1), y(2)})) ...
		; ...
            double(subs(cyclin_expr, {ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc}, ...
                {params.ksynth params.acdc25 params.bcdc25 params.ncdc25 params.awee1 params.bwee1 params.nwee1 params.adeg params.bdeg params.ndeg ...
                consts.EC50cdc25 consts.EC50wee1 consts.EC50deg y(1), y(2)}))
        ];
end


% function [cdk1, cyclin] = cdk1_cyclin_ode_nullclines(params, consts)
% 
%     % solve ckd1 rate eqn for nullcline.
%     syms ksynth;
%     cdk1_eqn = ksynth ...
% 	    + (acdc25 + params.bcdc25*(y(1)^params.ncdc25)/(consts.EC50cdc25^params.ncdc25 + y(1)^params.ncdc25))*(y(2)-y(1)) ...
% 	    - (params.awee1 + params.bwee1*(consts.EC50wee1^params.nwee1)/(consts.EC50wee1^params.nwee1 + y(1)^params.nwee1))*y(1) ...
% 	    - (params.adeg + params.bdeg*(y(1)^params.ndeg)/(consts.EC50deg^params.ndeg + y(1)^params.ndeg))*y(1);
% 
%     cdk1 = 
% end


function [cdk1_expr, cyclin_expr] = cdk1_cyclin_ode_exprs
    syms ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc;
   
    % sym_exprs is a vector of rate equations for [Cdk1; Cyclin].
    cdk1_expr = ksynth ...
		+ (acdc25 + bcdc25*(Cdk1^ncdc25)/(EC50cdc25^ncdc25 + Cdk1^ncdc25))*(Cyc-Cdk1) ...
		- (awee1 + bwee1*(EC50wee1^nwee1)/(EC50wee1^nwee1 + Cdk1^nwee1))*Cdk1 ...
		- (adeg + bdeg*(Cdk1^ndeg)/(EC50deg^ndeg + Cdk1^ndeg))*Cdk1;

    cyclin_expr = ksynth ...
        - (adeg + bdeg*(Cdk1^ndeg)/(EC50deg^ndeg + Cdk1^ndeg))*Cyc;
end