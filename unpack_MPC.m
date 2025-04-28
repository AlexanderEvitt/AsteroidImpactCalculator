function observations = unpack_MPC(filename,start_epoch)
    % Reads MPC observation data
    % https://www.minorplanetcenter.net/iau/info/OpticalObs.html
    % Read file
    lines = readlines(filename);

    % Disable time adjustment warning
    juliandate(2024,10,17.1);
    w = warning('query','last');
    id = w.identifier;
    warning('off',id);

    % Initialize observations table (time, ra,dec)
    % Time is seconds since JD TT = 2460600.5
    % UTC + 37s + 32.184s
    n = length(lines);
    observations = zeros(n,3);

    % Iterate and add to table
    for l = 1:n-1
        % Get specific line in text
        line = char(lines(l));

        % Process time
        year = str2double(line(16:19));
        month = str2double(line(21:22));
        day = str2double(line(24:32));
        JD = juliandate(year,month,day) - start_epoch; % time since epoch
        t = JD*86400 + 32.184 + 37; % seconds since 17 Oct 2024 (epoch from MPC)

        % Process right ascension
        ra.hour = str2double(line(33:34));
        ra.min = str2double(line(36:37));
        ra.sec = str2double(line(39:44));

        ra.deg = 15*ra.hour + 0.25*ra.min + (15/3600)*ra.sec;

        % Process declination
        dec.deg = str2double(line(45:47));
        dec.min = str2double(line(49:50));
        dec.sec = str2double(line(52:56));

        % Get sign of degrees
        if dec.deg == 0
            sign = 1;
        else
            sign = dec.deg/abs(dec.deg);
        end

        % Convert to degrees fully
        dec.deg = dec.deg + sign*(1/60)*dec.min + sign*(1/3600)*dec.sec;
        
        % Append to table
        observations(l,:) = [t,ra.deg,dec.deg];
    end
end