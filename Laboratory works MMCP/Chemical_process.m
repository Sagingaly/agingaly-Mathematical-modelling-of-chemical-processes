n = 10;
x_0 = 1/2 + 1/(n+3);
r = 1 + 2 / (n+3);
k = 1 + 1 / (n+1);

fprintf('x_0 value is: %f\n', x_0);
fprintf('r value is: %f\n', r);
fprintf('k value is: %f\n', k);

t1 = 1;
t2 = 2;

x1 = exp(r*t1) / ((1/x_0) + (1/k) * (exp(r*t1) - 1));
x2 = exp(r*t2) / ((1/x_0) + (1/k) * (exp(r*t2) - 1));

fprintf('x1 value is: %f\n', x1);
fprintf('x2 value is: %f\n', x2);

epsilon1 = 1e-3;
epsilon2 = 1e-4;

x1_adj = x1 * (1 + epsilon1);
x2_adj = x2 * (1 + epsilon2);

fprintf('x1_adj_is: %f\n', x1_adj);
fprintf('x2_adj_is: %f\n', x2_adj);

k_bar = ((exp(r * t1) / x1_adj) - (exp(r * t2) / x2_adj)) / (exp(r*t1) - exp(r*t2));  
x_0_bar = (exp(r * t1) / x1_adj) - k_bar * (exp(r * t1) - 1);

fprintf('k_bar_is: %f\n', k_bar);
fprintf('x_0_is: %f\n', x_0_bar);

k_bar_bar = 1 / k_bar;
x_0_bar_bar = 1 / x_0_bar;

delta_k = abs(k - k_bar_bar);
delta_x_0 = abs(x_0 - x_0_bar_bar);

fprintf('Delta_k: %f\n', delta_k);
fprintf('Delta_X: %f\n', delta_x_0);

% Создаем временной диапазон
t_values = linspace(0, 5, 100);

% Создаем массивы для delta_k и delta_x_0, содержащие одинаковые значения для всех t
delta_k_values = delta_k * ones(size(t_values));
delta_x_0_values = delta_x_0 * ones(size(t_values));

% Построение графика
figure;
plot(t_values, delta_k_values, 'b-', 'LineWidth', 2);
hold on;
plot(t_values, delta_x_0_values, 'r--', 'LineWidth', 2);
hold off;

title('Comparing analytical and theoretical values');
xlabel('Time');
ylabel('Values of function');
legend('Delta k', 'Delta x_0');
grid on;
t1 = 1;              
t2 = 
eps = 0.4;            

k_analytic = 1 / (3 * sqrt(n + 1));
A_analytic = n * (n + 1)^(1/3);