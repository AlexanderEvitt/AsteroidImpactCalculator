function dYdt = combined_dynamics(tc, Y, F, n)
    %% Combines the dynamics of F and Phi together
    % Allows integration of both together
    X = Y(1:n);                % nominal state
    Phi = reshape(Y(n+1:end), n, n); % reshape flattened Phi
    
    % Evaluate dynamics
    dXdt = F(tc, X);
    
    % Evaluate A matrix
    A = numerical_jacobian(F, X, tc);
    
    % Evaluate Phi_dot
    dPhidt = A * Phi;
    
    % Return flattened form
    dYdt = [dXdt; dPhidt(:)];
end