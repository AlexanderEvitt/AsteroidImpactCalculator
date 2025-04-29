function ephemeris = pack_ephemeris(JDs)
    % Check if ephemeris data is present in repo
    if isfile("ephemeris.mat")
        load("ephemeris.mat");
        ephemeris = ephemeris;
    else
        % Puts all the ephemeris data into a struct
        ephemeris.dates = JDs;
        ephemeris.sun = planetEphemeris(JDs,"SolarSystem","Sun","432t").';
        ephemeris.venus = planetEphemeris(JDs,"SolarSystem","Venus","432t").';
        ephemeris.earth = planetEphemeris(JDs,"SolarSystem","Earth","432t").';
        ephemeris.moon = planetEphemeris(JDs,"SolarSystem","Moon","432t").';
        ephemeris.mars = planetEphemeris(JDs,"SolarSystem","Mars","432t").';
        ephemeris.jupiter = planetEphemeris(JDs,"SolarSystem","Jupiter","432t").';
        % From https://ssd.jpl.nasa.gov/astro_par.html
        ephemeris.mu_sun = 1.32712440041279419e11;
        ephemeris.mu_venus = 324858.592000;
        ephemeris.mu_earth = 398600.435507;
        ephemeris.mu_moon = 4902.800118;
        ephemeris.mu_mars = 42828.375816; % system, not planet
        ephemeris.mu_jupiter = 126712764.100000; % system, not planet
    end
end