function dydt = solar_system_force_model(t,y,ephemeris,start_epoch)
    % Get Julian date (time defined in julian seconds since epoch)
    JD = start_epoch + (t/86400);
    
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
    %r = x - sun_pos;
    %er = r/norm(r);
    %ev = v/norm(v);
    %c = 299792.458;
    %a = a + (ephemeris.mu_sun/(norm(r)^2))*((4*(ephemeris.mu_sun/(norm(r)*c^2)) - ((norm(v)^2)/(c^2)))*er + 4*((norm(v)^2)/(c^2))*dot(er,ev)*ev);


    % Add Yarkovsky effect
    %a = a-(2e-12/1000)*v/norm(v); % empirical estimate

    % Pass to integrator
    dydt(1:3,1) = v;
    dydt(4:6,1) = a;
end

function acc = gravity(mu, r)
    r2 = dot(r, r);          % r squared (no sqrt)
    r3 = r2 * sqrt(r2);      % r cubed (1 sqrt + 1 multiply)
    acc = -(mu / r3) * r;
end