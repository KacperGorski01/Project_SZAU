%% 
clear; clc;

%% Porównanie modelu nieliniowego i zlinearizowanego

Tend = 20000;

% Nieliniowy model dynamiczny

function f = fNonlinear(t,h)
    tau = 125;
    F1 = 200 - 100*(t-tau >= 1000) + 100*(t-tau >= 10000);
    FD = 100 - 50*(t >= 6000) + 50*(t >= 10000);
    f1 = ( F1 + FD - 23 * sqrt(h(1))) / (0.7 * h(1));
    f2 = ( 23 * sqrt(h(1)) - 30 * sqrt(h(2))) / (1.35 * h(2)^2);
    f = [f1; f2];
end

h0 = [ 170.1323 ; 100 ];

[t,h] = ode45(@fNonlinear, [0, Tend], h0, odeset('RelTol',1e-6));

% Zlinearizowany model dynamiczny
% x1 = delta_h1 = h1 - 170.1323
% x2 = delta_h2 = h2 - 100
% u1 = delta_F1 = F1 - 200
% u2 = delta_FD = FD - 100

function f = fLinear(t,x)
    tau = 125;
    A = [-0.0074032, 0; 0.0000653, -0.00011111];
    B = [0.00839683, 0.00839683; 0, 0];
    u1 = -100*(t-tau >= 1000) + 100*(t-tau >= 10000);
    u2 = -50*(t >= 6000) + 50*(t >= 10000);
    f = A*x + B*[u1; u2];
end

x0 = [0; 0];
[t2,x] = ode45(@fLinear, [0, Tend], x0, odeset('RelTol',1e-6));
h1 = x(:,1) + 170.1323;
h2 = x(:,2) + 100;


figure(1);
sgtitle({
    '$F_{in} = 200 - 100 \cdot (t \ge 1000) + 100 \cdot (t \ge 10000)$', ...
    '$F_{D} = 100 - 50 \cdot (t \ge 6000) + 50 \cdot (t \ge 10000)$'
    }, 'Interpreter','latex','FontSize',14);
subplot(2,1,1);
hold on;
plot(t, h(:,1), 'b', 'LineWidth', 1.5); 
plot(t2, h1, 'r--', 'LineWidth', 1.5);
legend('Model nieliniowy', 'Model zlinearizowany', 'Location', 'southeast');
grid on
grid minor
subplot(2,1,2);
hold on;
plot(t, h(:,2), 'b', 'LineWidth', 1.5); 
plot(t2, h2, 'r--', 'LineWidth', 1.5);
legend('Model nieliniowy', 'Model zlinearizowany', 'Location', 'southeast');
grid on;
grid minor
xlabel('Czas [s]');
ylabel('Poziom cieczy [cm]');
