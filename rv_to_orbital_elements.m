function [a,e,i,RAAN,argp,theta] = rv_to_orbital_elements(r,v,mu)
    % Calc angular momentum vector
    h = cross(r,v);
    n = cross([0,0,1],h);

    % Calc eccentricity vector and mag
    e_bar = ((norm(v)^2-mu/norm(r))*r-dot(r,v)*v)/mu;
    e = norm(e_bar);

    specific_energy = norm(v)^2/2-mu/norm(r);
    
    % Presuming orbit is an ellipse
    a = -mu/(2*specific_energy);
    p = a*(1-e^2);

    i = acos(h(3)/norm(h));

    RAAN = acos(n(1)/norm(n));

    if n(2)<0
       RAAN = (2*pi)-RAAN;
    end
    
    argp = acos(dot(n,e_bar)/(norm(n)*e));
    
    if e_bar(3)<0
       argp = (2*pi)-argp;
    end

    theta = acos( dot(e_bar,r)/(e*norm(r)) );
    
    if dot(r,v)<0
       theta = (2*pi) - theta;
    end
    
    % Convert back to degrees
    i = (180/pi)*i;
    argp = (180/pi)*argp;
    RAAN = (180/pi)*RAAN;
    theta = (180/pi)*theta;
end