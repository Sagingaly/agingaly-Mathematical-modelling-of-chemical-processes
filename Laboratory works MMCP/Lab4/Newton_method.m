% Initialize parameters
x0 = 1/2 + 4/12;
r_initial = 0.3 + 1/17;  % Initial value of r
r = r_initial;
k = 1 + 1/17;


t1 = 1;
t2 = 3;
t3 = 5;

% Calculate x1, x2, x3 based on logistic growth model
x1 = exp(t1*r) / ((1/x0) + 1/k * (exp(t1*r) - 1));
x2 = exp(t2*r) / ((1/x0) + 1/k * (exp(t2*r) - 1));
x3 = exp(t3*r) / ((1/x0) + 1/k * (exp(t3*r) - 1));

% Display results
fprintf('x1 value is: %f\n', x1);
fprintf('x2 value is: %f\n', x2);
fprintf('x3 value is: %f\n', x3);

% Newton's method initialization
r0 = 0.358823;  % initial guess for Newton's method

% Define the function f(r) and its derivative f'(r)
f_r_sym = @(r_sym) (x2 * exp(r_sym * t1) - x1 * exp(r_sym * t2)) * (exp(r_sym * t3) - exp(r_sym * t1)) + ...
                   (exp(r_sym * t1) - exp(r_sym * t2)) * (x2 * exp(r_sym * t1) - (x1 * x2 / x3) * exp(r_sym * t3));

f_prime_r_sym = @(r_sym) (x2 * t1 * exp(r_sym * t1) - x1 * t2 * exp(r_sym * t2)) * (exp(r_sym * t3) - exp(r_sym * t1)) + ...
                        (x2 * exp(r_sym * t1) - x1 * exp(r_sym * t2)) * (t3 * exp(r_sym * t3) - t1 * exp(r_sym * t1)) + ...
                        (t1 * exp(r_sym * t1) - t2 * exp(r_sym * t2)) * (x2 * t1 * exp(r_sym * t1) - (x1 * x2 / x3) * t3 * exp(r_sym * t3));

% Newton's method iteration
max_iterations = 100;
tolerance = 1e-6;
for i = 1:max_iterations
    f_val = f_r_sym(r0);
    f_prime_val = f_prime_r_sym(r0);
    
    if abs(f_prime_val) < tolerance
        error('Derivative near zero, method fails to converge.');
    end
    
    % Newton's update
    r_new = r0 - f_val / f_prime_val;
    
    % Check convergence
    if abs(r_new - r0) < tolerance
        fprintf('Newton method converged after %d iterations\n', i);
        break;
    end
    
    r0 = r_new;  % update for next iteration
end

% Final value of r after convergence
fprintf('Final value of r after Newton method: %f\n', r0);


% Calculate error Delta_r = |r - r_n+1|
delta_r = abs(r_initial - r0);
fprintf('Error Delta_r: %f\n', delta_r);

% Introduce small error E for adjusted values
E = 1e-4;

x1_adj = x1 * (1 + E);
x2_adj = x2 * (1 + E);
x3_adj = x3 * (1 + E);

% Display adjusted values
fprintf('x1_adj value is: %f\n', x1_adj);
fprintf('x2_adj value is: %f\n', x2_adj);
fprintf('x3_adj value is: %f\n', x3_adj);

% Calculate k_bar_bar and x_0_bar_bar using adjusted values
k_bar_bar = ((exp(r * t1) / x1_adj) - (exp(r * t2) / x2_adj)) / (exp(r*t1) - exp(r*t2));  
x_0_bar_bar = (exp(r * t1) / x1_adj) - k_bar_bar * (exp(r * t1) - 1);

% Correct the calculation of k_bar and x_0_bar
k_bar = 1 / k_bar_bar;
x_0_bar = 1 / x_0_bar_bar;

% Calculate the differences (errors) between original and adjusted parameters
delta_k = abs(k - k_bar);
delta_x_0 = abs(x0 - x_0_bar);

% Display delta values
fprintf('Delta_k: %f\n', delta_k);
fprintf('Delta_x_0: %f\n', delta_x_0);

% Create a time range for plotting
t_values = linspace(0, 5, 100);

% Create arrays for delta_k and delta_x_0 with constant values for all t
delta_k_values = delta_k * ones(size(t_values));
delta_x_0_values = delta_x_0 * ones(size(t_values));

% Plotting
figure;
plot(t_values, delta_k_values, 'b-', 'LineWidth', 2);
hold on;
plot(t_values, delta_x_0_values, 'r--', 'LineWidth', 2);
hold off;

title('Comparison of analytical and theoretical values');
xlabel('Time');
ylabel('Values');
legend('Delta k', 'Delta x_0');
grid on;
