function [r, v] = orbital_elements_to_rv(a,e,i,RAAN,argp,theta,mu)
    % Convert angular parameters from degrees into radians
    i  = (pi/180)*i;
    RAAN = (pi/180)*RAAN;
    argp = (pi/180)*argp;
    theta = (pi/180)*theta;

    % Calculate r, vr, vt using the orbit parameter p
    p = a*(1-e^2);
    r = p/(1+e*cos(theta));
    vr = sqrt(mu/p)*e*sin(theta);
    vt = sqrt(mu/p)*(1+e*cos(theta));

    % Do the necessary 3-1-3 rotation to orient everything right
    % The matrix transformations are slightly different in convention because I'm
    % referencing some slides from another class lol
    Rom = [cos(RAAN) , -sin(RAAN) , 0; sin(RAAN), cos(RAAN), 0; 0, 0, 1];
    Ri = [1, 0, 0; 0, cos(i), -sin(i); 0, sin(i), cos(i)];
    Rthw = [cos(theta+argp) , -sin(theta+argp) , 0; sin(theta+argp), cos(theta+argp), 0; 0, 0, 1];

    % Construct the output
    r = Rom*Ri*Rthw*[r;0;0];
    v = Rom*Ri*Rthw*[vr;vt;0];
end