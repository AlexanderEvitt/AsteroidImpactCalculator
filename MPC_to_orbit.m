function X = MPC_to_orbit(a,e,i,RAAN,argp,M,n)
    % Returns the state of the body from the MPC ephemeris
    mu = 1.327124400189e11;
    M = deg2rad(M);
    
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