function J = numerical_jacobian(func, X, t, epsilon)
    %% Numerical Jacobian algorithm
    % Written by Alexander Evitt, 5 Feb 2025
    % func: Function handle
    % X: State vector where Jacobian is evaluated
    % epsilon: Perturbation size
    
    if nargin < 4
        epsilon = 1e-6;
    end
    
    n = length(X);
    f0 = func(t,X);
    m = length(f0);
    J = zeros(m, n);

    for i = 1:n
        Xperturb = X;
        Xperturb(i) = Xperturb(i) + epsilon; % Perturb one variable
        fPerturb = func(t,Xperturb);
        J(:, i) = (fPerturb - f0) / epsilon; % Finite difference approximation
    end
end