function pos = unpack_ephemeris(dates, ephemeris, JD)
    % Returns position of body given date and ephemeris data
    % Vectorized spline interpolation for all three components at once
    pos = interp1(dates.', ephemeris.', JD, 'spline').';
end