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

% Exposer function.
function funcs = cdk1_cyclin_ode
    funcs.dydt = @cdk1_cyclin_ode_dydt;
    funcs.nullclines = @cdk1_cyclin_ode_nullclines;
    funcs.dydt_array = @cdk1_cyclin_ode_dydt_array;
end


function dydt = cdk1_cyclin_ode_dydt(t, y, params, consts)

    syms ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc pos_fb_strength;

    % dydt is a vector of rates for [Cdk1; Cyclin].
    [cdk1_expr, cyclin_expr] = cdk1_cyclin_ode_exprs;
    dydt = [ double(subs(cdk1_expr, {ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc pos_fb_strength}, ...
                {params.ksynth params.acdc25 params.bcdc25 params.ncdc25 params.awee1 params.bwee1 params.nwee1 params.adeg params.bdeg params.ndeg ...
                consts.EC50cdc25 consts.EC50wee1 consts.EC50deg y(1) y(2) params.pos_fb_strength})) ...
		; ...
            double(subs(cyclin_expr, {ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc pos_fb_strength}, ...
                {params.ksynth params.acdc25 params.bcdc25 params.ncdc25 params.awee1 params.bwee1 params.nwee1 params.adeg params.bdeg params.ndeg ...
                consts.EC50cdc25 consts.EC50wee1 consts.EC50deg y(1) y(2) params.pos_fb_strength}))
        ];
end

function dydt = cdk1_cyclin_ode_dydt_array(t, y, params, consts)

    % serves the same purpose as `cdk1_cyclin_ode_dydt` but takes y as a 3D
    % array and returns dydt as a 3D array --> helper for generating vector field.
    syms ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc pos_fb_strength;

    % dydt is a vector of rates for [Cdk1; Cyclin].
    [cdk1_expr, cyclin_expr] = cdk1_cyclin_ode_exprs;
    dydt = double(subs(cdk1_expr, {ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc pos_fb_strength}, ...
                {params.ksynth params.acdc25 params.bcdc25 params.ncdc25 params.awee1 params.bwee1 params.nwee1 params.adeg params.bdeg params.ndeg ...
                consts.EC50cdc25 consts.EC50wee1 consts.EC50deg y(:, :, 1) y(:, :, 2) params.pos_fb_strength}));
    % append to third axis.
	dydt = cat(3, dydt, double(subs(cyclin_expr, {ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc pos_fb_strength}, ...
                {params.ksynth params.acdc25 params.bcdc25 params.ncdc25 params.awee1 params.bwee1 params.nwee1 params.adeg params.bdeg params.ndeg ...
                consts.EC50cdc25 consts.EC50wee1 consts.EC50deg y(:, :, 1) y(:, :, 2) params.pos_fb_strength})));
end


function [cdk1_val_pts, cyclin_val] = cdk1_cyclin_ode_nullclines(params, consts, cyc_vals, hyst_range)

	syms ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc pos_fb_strength;
	syms rate;

	% set up helpers.
	null_eqn = rate == 0;
	[cdk1_expr, cyclin_expr] = cdk1_cyclin_ode_exprs;
    
	% solve rate eqns for nullclines.
	% - Cdk1.
    % numerically solve for each Cyc value in cyc_vals.
    cdk1_val_pts = zeros([length(cyc_vals) 2]);
    for k=1:length(cyc_vals)
        cdk1_eqn = compose(null_eqn, cdk1_expr);
        cdk1_eqn = subs(cdk1_eqn, {ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cyc pos_fb_strength}, ...
                    {params.ksynth params.acdc25 params.bcdc25 params.ncdc25 params.awee1 params.bwee1 params.nwee1 params.adeg params.bdeg params.ndeg ...
                    consts.EC50cdc25 consts.EC50wee1 consts.EC50deg cyc_vals(k) params.pos_fb_strength});
        cdk1_val_pts(k, 2) = vpasolve(cdk1_eqn, Cdk1);
        cdk1_val_pts(k, 1) = cyc_vals(k);
    end

    % increment intervals for searching multiple solutions for Cdk1 in
    % hysteretic region -- smaller value => multiple solns very close
    % together are missed.
    hyst_dx = 2;
    % attempt to find multiple solns if Cyc lies in hyst_range.
    addnl_k = k+1;
    for addnl_val=hyst_range(1):1:hyst_range(2)

        % set up the equation.
        cdk1_eqn = compose(null_eqn, cdk1_expr);
        cdk1_eqn = subs(cdk1_eqn, {ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cyc pos_fb_strength}, ...
                    {params.ksynth params.acdc25 params.bcdc25 params.ncdc25 params.awee1 params.bwee1 params.nwee1 params.adeg params.bdeg params.ndeg ...
                    consts.EC50cdc25 consts.EC50wee1 consts.EC50deg addnl_val params.pos_fb_strength});
        
        % get first soln; closest to zero.
        first_soln = vpasolve(cdk1_eqn, Cdk1, 0);
        
        % now move soln intervals by hyst_dx to search.
        solns_cnt = 1;
        for soln_start=first_soln+hyst_dx:hyst_dx:max(cyc_vals)
            soln = vpasolve(cdk1_eqn, Cdk1, [soln_start soln_start+hyst_dx]);
            if not(isempty(soln))
                cdk1_val_pts(addnl_k, 2) = soln;
                cdk1_val_pts(addnl_k, 1) = addnl_val;
                solns_cnt = solns_cnt + 1;
                addnl_k = addnl_k + 1;
            end
            if solns_cnt>2
                break;
            end
        end
    end
    
	% - Cyclin.
	cyclin_eqn = compose(null_eqn, cyclin_expr);
    cyclin_eqn = subs(cyclin_eqn, {ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg}, ...
                {params.ksynth params.acdc25 params.bcdc25 params.ncdc25 params.awee1 params.bwee1 params.nwee1 params.adeg params.bdeg params.ndeg ...
                consts.EC50cdc25 consts.EC50wee1 consts.EC50deg});
	cyclin_val = solve(cyclin_eqn, Cyc);

end


% Helper function that renders the ODEs as symbolic equations.
function [cdk1_expr, cyclin_expr] = cdk1_cyclin_ode_exprs
    syms ksynth acdc25 bcdc25 ncdc25 awee1 bwee1 nwee1 adeg bdeg ndeg EC50cdc25 EC50wee1 EC50deg Cdk1 Cyc pos_fb_strength;
   
    % sym_exprs is a vector of rate equations for [Cdk1; Cyclin].
    cdk1_expr = ksynth ...
		+ pos_fb_strength*(acdc25 + bcdc25*(Cdk1^ncdc25)/(EC50cdc25^ncdc25 + Cdk1^ncdc25))*(Cyc-Cdk1) ...
		- pos_fb_strength*(awee1 + bwee1*(EC50wee1^nwee1)/(EC50wee1^nwee1 + Cdk1^nwee1))*Cdk1 ...
		- (adeg + bdeg*(Cdk1^ndeg)/(EC50deg^ndeg + Cdk1^ndeg))*Cdk1;

    cyclin_expr = ksynth ...
        - (adeg + bdeg*(Cdk1^ndeg)/(EC50deg^ndeg + Cdk1^ndeg))*Cyc;
end