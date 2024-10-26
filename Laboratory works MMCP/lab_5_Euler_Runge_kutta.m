N = 10;           
tmax = 5;       
r = 0.002;        
alfa = 0.003;    
b = 0.0021;       
tau = 0.001;      
x_i = 0.0032;

t_delta = tmax/N;
t_values = linspace(0, tmax, N+1);


eu_values = zeros(1, N+1);
rk_values = zeros(1, N+1);

eu_values(1) = x_i;
rk_values(1) =x_i;

for i = 1:N
    x_derivative = -r*x_i + (alfa * b * x_i^2) / b + tau * x_i;
   
    x_i1 = x_i + t_delta* 0.5 * x_derivative;
    
    rk_values(i+1) = x_i1 

    x_i = x_i1

end


x_i = 0.0032;


for i = 1:N
    x_derivative = -r*x_i + (alfa * b * x_i^2) / b + tau * x_i;
   
    x_i_1 = x_i + t_delta * x_derivative;
    
    eu_values(i+1) = x_i_1 

    x_i = x_i_1

end

figure;
plot(t_values, eu_values, 'r-', LineWidth= 2);
hold on;
plot(t_values, rk_values, 'b--', LineWidth= 2);
hold off;


title('Comparison between Euler and Runge- Kutta method');
xlabel('Time');
ylabel('Values');
legend('Euler','Runge- Kutta');
grid on;