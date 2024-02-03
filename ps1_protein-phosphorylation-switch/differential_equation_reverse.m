function dxdt = differential_equation_reverse(t, x, K, P, xM, n)
    % Return the rate of change of x as a function of x and t
    % and the parameters of the equation.
    dxdt = P*x - (K.*(1 - x).*(((x.^n)./(x.^n + xM^n))+0.1));
end