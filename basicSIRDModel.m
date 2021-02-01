 

% System of ODEs
f = @(i, s) [-0.9*s(1)*s(2); 0.9*s(1)*s(2)-0.3*s(2)-0.005*s(2); 0.3*s(2); 0.005*s(2)];
 
% Intervals
a = 0;
b = 150;

% Step size (in days)
h = 1;
 
% Initial conditions
ya = [0.999; 0.001; 0; 0];

[t, y] = rungaKuttaSystem(f, a, b, ya, h);
 
figure
plot(t, y);
title("SIRD Model of Infectious Disease Spread")
legend("Suspectible", "Infected", "Recovered", "Deceased");
xlabel("Days")
ylabel("Percentage")