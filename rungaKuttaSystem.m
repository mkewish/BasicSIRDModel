% ==================================================
% 4th Order Runga-Kutta Method for Solving Systems
% ==================================================
% Inputs:
% f = function handle corres. to ODE y' = f(t,y) 
% a = starting t-value for solution interval
% b = max t-value for solution interval 
% ya = vector of y(a) (initial condition vector)
% h = step size
% 
% Outputs:
% t = vector of t-values at which soln is computed
% y = approximated soln. at each t-value

function [t, y] = rungaKuttaSystem(f, a, b, ya, h)
    
    neq = numel(ya);
    
    N = floor((b-a)/h);
    t = zeros(N+1, 1)';
    y = zeros(N+1, neq)'; 
    
    t(1) = a;
    y(:, 1) = ya;
    
    for ii = 2:N+1                              
        k_1 = f(t(ii-1), y(:, ii-1));
        k_2 = f(t(ii-1)+0.5*h, y(:, ii-1)+(0.5*h*k_1));
        k_3 = f((t(ii-1)+0.5*h), (y(:, ii-1)+0.5*h*k_2));
        k_4 = f((t(ii-1)+h), (y(:, ii-1)+k_3*h));
        y(:, ii) = y(:, ii-1) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
        t(ii) = t(ii-1) + h;
    end
    
end