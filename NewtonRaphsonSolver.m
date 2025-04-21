function sol = NewtonRaphsonSolver(init, tol, f, f_prime)
    % Create placeholders for loop
    diff = 1;
    prev_sol = init;
    i = init;
    f_at_sol = f(init);
    
    % Iterate until tolerance is met
    while(diff >= tol) 
        % Apply Newton method
       sol = prev_sol - f(prev_sol)/f_prime(prev_sol);

       % Update values
       diff =  abs(sol - prev_sol);
       f_at_sol = [f_at_sol, f(prev_sol)];
       i = [i, sol];
       prev_sol = sol;
    end
end