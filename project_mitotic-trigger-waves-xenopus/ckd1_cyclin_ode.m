%{
Parameter	Value
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

function dydt = ckd1_cyclin_ode(t, y, params, consts)

    % dydt is a vector for [Cdk1; Cyclin].

    dydt = [ ...
            params.ksynth ...
			+ (params.acdc25 + params.bcdc25*(y(1)^params.ncdc25)/(consts.EC50cdc25^params.ncdc25 + y(1)^ncdc25))*(y(2)-y(1)) ...
			- (params.awee1 + params.bwee1*(consts.EC50wee1^params.nwee1)/(consts.EC50wee1^params.nwee1 + y(1)^params.nwee1))*y(1) ...
			- (params.adeg + params.bdeg*(y(1)^params.ndeg)/(consts.EC50deg^params.ndeg + y(1)^params.ndeg))*y(1) ...
		; ...
            params.ksynth ...
            - (params.adeg + params.bdeg*(y(1)^params.ndeg)/(consts.EC50deg^params.ndeg + y(1)^params.ndeg))*y(2)
        ]

end