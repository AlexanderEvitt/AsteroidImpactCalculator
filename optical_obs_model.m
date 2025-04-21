function y = optical_obs_model(t,X,ephemeris)
    % Measurement model, converts time and state to azi,elev
    c = 299792.458;

    % Get Julian date (time defined in julian seconds since 2000)
    JD = juliandate(2000,1,1,0,0,0) + (t/86400);

    % RA and DEC are independent of observer position on Earth
    % Also independent of Earth rotation
    % All that matters is where Earth is (R)
    R = unpack_ephemeris(ephemeris.dates,ephemeris.sun,JD);

    tau = 0;
    n = 5;
    for i = 1:n
        % Get body position using Keplerian approximation
        Xt = central_force_motion(1.327124400189e11,X,-tau);
        r = Xt(1:3);
    
        % Iterate for light-time (assumed downlink)
        tau = (1/c)*norm(r - R);
    end

    % Model position of body at t - tau
    Xt = central_force_motion(1.327124400189e11,X,-tau);
    r = Xt(1:3);

    % Convert barycentric to Earth-centered
    r_ECI = r - R;
    r_hat = r_ECI/norm(r_ECI);

    % Declination in radians
    dec = asin(r_hat(3));  % z component gives declination
    
    % Right ascension in radians
    ra = atan2(r_hat(2), r_hat(1));
    if ra < 0
        ra = ra + 2*pi;
    end

    % Calculate RA/DEC
    ra = rad2deg(ra);
    dec = rad2deg(dec);
    y = [ra; dec];
end


function pos = unpack_ephemeris(dates,ephemeris,JD)
    % Returns position of body given date and ephemeris data
    pos1 = interpn(dates.',ephemeris(1,:),JD, "spline");
    pos2 = interpn(dates.',ephemeris(2,:),JD, "spline");
    pos3 = interpn(dates.',ephemeris(3,:),JD, "spline");
    pos = [pos1;pos2;pos3];
end