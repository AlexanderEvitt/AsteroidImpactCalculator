function [X0, nominal, P0] = batch_estimate(F,G,t,Y,P,R,iters,X0,offset)
    %% Batch estimation algorithm
    % Written by Alexander Evitt, 7 Feb 2025
    % F: dynamics function handle
    % G: measurement model function handle
    % t: m x 1 matrix of times of measurements
    % Y: m x 1 matrix of values of measurements
    % R: n x n covariance matrix, W = inv(R)
    % X0: nominal initial state
    % offset : offset of 1 because there's no data at zero
    % iters: number of iterations

    % Set constants
    m = size(Y,1);
    n = size(X0,1);
    x0 = zeros(n,1);
    apriori = true;
    if nargin < 8 | apriori == false
        apriori = false;
        X0 = 100000*ones(n,1);
    end

    % Run iters number of times
    for k = 1:iters

        % Set from a priori
        if apriori
            Lambda = inv(P);
            N = inv(P)*x0;
        else
            Lambda = zeros(n);
            N = zeros(n,1);
        end

        % Integrate nominal trajectory
        options = odeset("RelTol",1e-12,"AbsTol",1e-14);
        [~,sols] = ode45(F,t,X0,options);

        % Iterate through measurements
        for i = 1+offset:m
            % Integrate nominal trajectory
            X_nom = sols(i,:).';

            % Recalculate A
            A = numerical_jacobian(F,X_nom,t(i));

            % Recalculate phi
            Phi = expm(A*(t(i) - t(1)));

            % Accumulate observation
            squigglyH = numerical_jacobian(G,X_nom,t(i));
            y = Y(i-offset,:).' - G(t(i),X_nom);
            H = squigglyH*Phi;
            Lambda = Lambda + (H.')*inv(R)*H;
            N = N + (H.')*inv(R)*y;
        end

        % Solve normal equations
        xhat = inv(Lambda)*N;

        % Update nominal trajectory
        X0 = X0 + xhat;

        % Shift
        apriori = true;
        x0 = x0 - xhat;
    end

    P0 = inv(Lambda);
    options = odeset("RelTol",1e-12,"AbsTol",1e-14);
    [~,sols] = ode45(F,t,X0,options);
    nominal = sols.';
end