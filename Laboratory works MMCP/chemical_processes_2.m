% Given values
E = 10e-4;  % Error term
n = 30;  % Number of time values
A = 1;  % Initial A value (adjust as needed)
k = 0.5;  % Initial k value (adjust as needed)

% Time values (t1, t2, ..., t30)
t = linspace(1, 30, n);  % Example time values from t1 to t30

% Step 1: Calculate x_i for each t_i
x = A * exp(k * t);

% Step 2: Calculate adjusted values x_i_adj with error
x_adj = x .* (1 + E);

% Step 3: Compute required sums
sum_ti = sum(t);
sum_ti2 = sum(t.^2);
sum_ln_xi = sum(log(x));
sum_ln_xi_ti = sum(log(x) .* t);

% Step 4: Compute determinants for Delta, Delta_k, and Delta_lnA

% Delta matrix: [ sum(ti), sum(ti^2); n, sum(ti) ]
Delta = (sum_ti * sum_ti) - (n * sum_ti2);

% Delta_k matrix: [ sum(ln(xi) * ti), sum(ti^2); sum(ln(xi)), sum(ti) ]
Delta_k = (sum_ln_xi_ti * sum_ti) - (sum_ln_xi * sum_ti2);

% Delta_lnA matrix: [ sum(ln(xi) * ti), sum(ti); sum(ln(xi)), n ]
Delta_lnA = (sum_ln_xi_ti * n) - (sum_ln_xi * sum_ti);

% Step 5: Find k_bar and ln(A_bar) using the formulas
k_bar = Delta_k / Delta;
ln_A_bar = Delta_lnA / Delta;

% Step 6: Find the adjusted A_bar
A_bar = exp(ln_A_bar);

% Step 7: Find Delta values
deltaaa_1 = abs(k - k_bar);  % Delta_1 = abs(k - k_bar)
deltaaa_2 = abs(A - A_bar);  % Delta_2 = abs(A - A_bar)

% Display results
fprintf('Delta: %f\n', Delta);
fprintf('Delta_k: %f\n', Delta_k);
fprintf('Delta_lnA: %f\n', Delta_lnA);
fprintf('k_bar: %f\n', k_bar);
fprintf('ln(A_bar): %f\n', ln_A_bar);
fprintf('A_bar: %f\n', A_bar);
fprintf('Delta_1 (k difference): %f\n', deltaaa_1);
fprintf('Delta_2 (A difference): %f\n', deltaaa_2);

% Optional: Plot x values and adjusted x values
figure;
plot(t, x, 'b-', 'LineWidth', 2);
hold on;
plot(t, x_adj, 'r--', 'LineWidth', 2);
title('Original and Adjusted x(t) values');
xlabel('Time (t)');
ylabel('x(t)');
legend('Original x(t)', 'Adjusted x_{adj}(t)');
grid on;
