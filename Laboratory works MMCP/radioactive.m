%{
n = 10;
eps = 10e-4;
fprintf('eps=%.4e\n',eps);

k = 1/3*sqrt(n+1);
A = n * (n+1)^1/3;

t1 = 1;
t2 = 4;

function result = func_with_E(t, A, k, eps);
 result = A * exp(k*t)*(1+eps);
end 


x1=func_with_E(t1, A, k, eps);
x2=func_with_E(t2, A, k, eps);

K=log(x1/x2)/(t1-t2);
A1=x1/exp(K*t1);
A2=x2/exp(K*t2);
fprintf('Analytical k= %.4f,A = %.4f\n', k, A);
fprintf('Theoretical k = %.4f, A= %.4f\n', k, A1);
fprintf('Accuracy E_k = %.4f',abs(k-k), abs(A-A1));
%}


%{

n = 10;
eps = 10e-4;
fprintf('eps = %.4e\n', eps);

% Вычисление констант k и A
k = 1 / (3 * sqrt(n + 1));
A = n * (n + 1)^(1/3);

% Временные точки
t1 = 1;
t2 = 4;

% Определение функции
function result = func_with_E(t, A, k, eps)
    result = A * exp(k * t) * (1 + eps);
end

% Вычисление значений функции в точках t1 и t2
x1 = func_with_E(t1, A, k, eps);
x2 = func_with_E(t2, A, k, eps);

% Вычисление теоретических значений
K = log(x1 / x2) / (t1 - t2);
A1 = x1 / exp(K * t1);
A2 = x2 / exp(K * t2);

% Вывод аналитических и теоретических значений
fprintf('Analytical k = %.4f, A = %.4f\n', k, A);
fprintf('Theoretical k = %.4f, A = %.4f\n', K, A1);

% Вычисление и вывод точности
fprintf('Accuracy E_k = %.4f, E_A = %.4f\n', abs(k - K), abs(A - A1));
%}


% Параметры задачи
n = 10;               % Количество шагов
t1 = 1;               % Первый временной момент
t2 = 4;               % Второй временной момент
eps = 0.4;            % Погрешность

% 1) Вычисление аналитического k и A
k_analytic = 1 / (3 * sqrt(n + 1));
A_analytic = n * (n + 1)^(1/3);

% Вычисление X1 и X2 с eps
X1 = A_analytic * exp(k_analytic * t1) + eps;  % x1 + eps
X2 = A_analytic * exp(k_analytic * t2) + eps;  % x2 + eps

% 2) Теоретические значения k и A
k_theoretical = log(X1 / X2) / (t1 - t2);
A_theoretical = X1 / exp(k_theoretical * t1);

% 3) Вычисление погрешностей
delta_k = abs(k_theoretical - k_analytic);
delta_A = abs(A_theoretical - A_analytic);

% Вывод аналитических и теоретических значений
fprintf('Аналитическое k = %.10f, A = %.10f\n', k_analytic, A_analytic);
fprintf('Теоретическое k = %.10f, A = %.10f\n', k_theoretical, A_theoretical);
fprintf('Погрешность Δk = %.10f, ΔA = %.10f\n', delta_k, delta_A);

% 4) Построение графика аналитических и теоретических значений
t_values = linspace(0, 5, 100);  % Временные значения

% Функции для вычисления аналитических и теоретических значений
analytical_func = @(t) A_analytic * exp(k_analytic * t);  % Без eps
theoretical_func = @(t) A_theoretical * exp(k_theoretical * t);  % Без eps

% Значения функций
analytical_values = arrayfun(analytical_func, t_values);
theoretical_values = arrayfun(theoretical_func, t_values);

% Построение графика
figure;
plot(t_values, analytical_values, 'b-', 'LineWidth', 2);  % Синим - аналитическое
hold on;
plot(t_values, theoretical_values, 'r--', 'LineWidth', 2);  % Красным - теоретическое
hold off;

% Оформление графика
title('Сравнение аналитических и теоретических значений');
xlabel('Время');
ylabel('Значение функции');
legend('Аналитическое', 'Теоретическое');
grid on;
