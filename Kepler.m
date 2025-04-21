function X = Kepler(mu,X0,t)
    % Returns the position of a spacecraft following a Keplerian orbit
    % [t] seconds after the state X0 is provided

    % Convert r0,v0 into orbital elements
    r0 = X0(1:3);
    v0 = X0(4:6);
    [a,e,i,RAAN,argp,theta0] = rv_to_orbital_elements(r0,v0,mu)
    theta0 = deg2rad(theta0);

    % Returns position and velocity at [t] seconds past theta0
    E0 = 2*atan(tan(theta0/2)*sqrt((1-e)/(1+e)))
    M0 = E0 - e*sin(E0);
    
    n = sqrt(mu/a^3);
    M = M0 + n*t; % t seconds
    
    % Function handles for Kepler's equation
    f = @(E) E - e*sin(E) - M;
    f_prime = @(E) 1 - e*cos(E);
    
    % Apply Newton-Raphson method
    init = 1.2;
    tol = 1e-6;
    E = NewtonRaphsonSolver(init, tol, f, f_prime);
    theta = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

    theta = rad2deg(theta);
    
    % Get r, v
    [r, v] = orbital_elements_to_rv(a,e,i,RAAN,argp,theta,mu);
    X = [r;v];
end