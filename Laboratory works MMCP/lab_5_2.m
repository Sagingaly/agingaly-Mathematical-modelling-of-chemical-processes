% Runge- Kutta method of 2nd order and Euler method
 

% Given parameters
N = 10;           
tmax = 5;       
r = 0.002;        
alfa = 0.003;    
b = 0.0021;       
tau = 0.001;      
x_i = 0.0032;

t_delta = tmax/ N;

t_values = linspace(0, tmax, N+1);

for i = 1: N
    x_derivative = -r*x_i + (alfa * b * x_i^2) / b + tau * x_i;
   
    x_i1 = x_i + t_delta * 0.5 * x_derivative;

    disp(['x step at', num2str(i), 'is', num2str(x_i1)]);

    x_i = x_i1;
end

for i = 1: N
    x_derivative = -r*x_i + (alfa * b * x_i^2) / b + tau * x_i;
   
    x_i_1 = x_i + t_delta * x_derivative;

    disp(['x step at', num2str(i), 'is \n', num2str(x_i_1)]);

    x_i = x_i_1;
end



