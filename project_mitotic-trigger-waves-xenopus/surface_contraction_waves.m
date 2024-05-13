function funcs = surface_contraction_waves
    funcs.surface_intersection_eqns = @surface_intersection_eqns;
end


function [sphere_surface_eqn, wave_eqn] = surface_intersection_eqns
    syms x y v rc t x0 y0 t0;
    sphere_surface_eqn = y^2 + x^2 == rc^2;
    wave_eqn = (y-y0)^2 + (x-x0)^2 == v^2*(t-t0)^2;
end

