x = linspace(0, 10, 100);
y = trapmf(x, [2 4 6 8]);

plot(x, y);
xlabel('x');
ylabel('Membership Degree');
title('Trapezoidal Membership Function');
grid on
grid minor