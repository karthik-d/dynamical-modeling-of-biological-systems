function dxdt = differential_equation(t, x, K, P, xM, n)
    % Return the rate of change of x as a function of x and t
    % and the parameters of the equation.
    dxdt = K.*(1 - x).*(((x.^n)./(x.^n + xM^n))+0.1) - P*x;
end