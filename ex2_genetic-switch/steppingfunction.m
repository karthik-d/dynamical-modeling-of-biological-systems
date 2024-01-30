function r = steppingfunction(t)

% nasty, hard-coded function with no parameters

if t<1000
    r = floor(t/200);
elseif t<2000
    r = (5-floor((t-1000)/200));
else
    r = 0;
end
