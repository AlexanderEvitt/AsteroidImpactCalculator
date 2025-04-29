function X = central_force_motion(mu,X0,t)
    r = X0(1:3,:);
    v = X0(4:6,:);

    % Assume no forces during time
    r = r + v*t;
    X = [r;v];
end