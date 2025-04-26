function dydt = solar_system_force_model(t,y,ephemeris)
    % Get Julian date (time defined in julian seconds since epoch)
    JD = juliandate(2024,10,17,0,0,0) + (t/86400);
    
    % Get positions of Venus, Earth, Mars, Jupiter, and Sun
    sun_pos = unpack_ephemeris(ephemeris.dates,ephemeris.sun,JD);
    venus_pos = unpack_ephemeris(ephemeris.dates,ephemeris.venus,JD);
    earth_pos = unpack_ephemeris(ephemeris.dates,ephemeris.earth,JD);
    moon_pos = unpack_ephemeris(ephemeris.dates,ephemeris.moon,JD);
    mars_pos = unpack_ephemeris(ephemeris.dates,ephemeris.mars,JD);
    jupiter_pos = unpack_ephemeris(ephemeris.dates,ephemeris.jupiter,JD);

    % Unpack x and v, initialize a
    x = y(1:3);
    v = y(4:6);
    a = [0;0;0];

    % Add gravitational acceleration from each body
    a = a + gravity(ephemeris.mu_sun,x-sun_pos);
    a = a + gravity(ephemeris.mu_venus,x-venus_pos);
    a = a + gravity(ephemeris.mu_earth,x-earth_pos);
    a = a + gravity(ephemeris.mu_moon,x-moon_pos);
    a = a + gravity(ephemeris.mu_mars,x-mars_pos);
    a = a + gravity(ephemeris.mu_jupiter,x-jupiter_pos);

    % Add relativistic corrections for each body
    r = x - sun_pos;
    er = r/norm(r);
    ev = v/norm(v);
    c = 299792.458;
    %a = a + (ephemeris.mu_sun/(norm(r)^2))*((4*(ephemeris.mu_sun/(norm(r)*c^2)) - ((norm(v)^2)/(c^2)))*er + 4*((norm(v)^2)/(c^2))*dot(er,ev)*ev);


    % Add Yarkovsky effect
    %a = a-(2e-12/1000)*v/norm(v); % empirical estimate

    % Pass to integrator
    dydt(1:3,1) = v;
    dydt(4:6,1) = a;
end

function acc = gravity(mu,r)
    % Returns gravity in km/s^2
    acc = (-mu/(norm(r)^3))*r;
end

function pos = unpack_ephemeris(dates,ephemeris,JD)
    % Returns position of body given date and ephemeris data
    pos1 = interpn(dates.',ephemeris(1,:),JD, "spline");
    pos2 = interpn(dates.',ephemeris(2,:),JD, "spline");
    pos3 = interpn(dates.',ephemeris(3,:),JD, "spline");
    pos = [pos1;pos2;pos3];
end