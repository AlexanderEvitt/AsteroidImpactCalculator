function dydt = solar_system_force_model(t,y,ephemeris)
    % Get Julian date (time defined in julian seconds since 2000)
    JD = juliandate(2024,10,17,0,0,0) + (t/86400);
    
    % Get positions of Venus, Earth, Mars, Jupiter, and Sun
    sun_pos = unpack_ephemeris(ephemeris.dates,ephemeris.sun,JD);
    earth_pos = unpack_ephemeris(ephemeris.dates,ephemeris.sun,JD);

    % Unpack x and v, initialize a
    x = y(1:3);
    v = y(4:6);
    a = [0;0;0];

    % Add gravitational acceleration
    a = a + gravity(1.32724400e11,x-sun_pos);
    a = a + gravity(3.986004418e5,x-earth_pos);

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