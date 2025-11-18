%% Porównanie modelu nieliniowego i modelu SL
clear; clc;

Tend = 2e6;

h0 = [ 170 ; 100 ];
F1 = @(t) 200 + 100*(t>= 2e5) - 100*(t >= 1e6);
FD = @(t) 100 + 50*(t >= 5e5) - 50*(t >= 1.5e6);


%% ---------------------------------------------------------------
% Model nieliniowy

fNonlinear = @(t,h) [ ...
    ( F1(t) + FD(t) - 23 * sqrt(h(1))) / (0.7 * h(1)) ; ...
    ( 23 * sqrt(h(1)) - 30 * sqrt(h(2))) / (1.35 * h(2)^2) ...
];

[tn,hn] = ode45(fNonlinear, [0, Tend], h0, odeset('RelTol',1e-6));


%% ---------------------------------------------------------------
% Model SL

Ts = Tend/1e4;
tl = 0:Ts:Tend;

hl = zeros(length(tl), 2);
hl(1, :) = h0';

A = zeros(2,2);
B = zeros(2,2);
for i = 1:length(tl)-1
    t = tl(i);

    F1_p = F1(t);
    FD_p = FD(t);
    h1_p = ( 1/23 * (F1_p + FD_p) )^2;
    h2_p = ( 1/30 * (F1_p + FD_p) )^2;

    A(1,1) = -1.42857142857143*(F1_p + FD_p - 23*sqrt(h1_p))/h1_p^2 - 16.4285714285714/h1_p^(3/2);
    A(1,2) = 0;
    A(2,1) = 8.51851851851852/(sqrt(h1_p)*h2_p^2);
    A(2,2) = -1.48148148148148*(23*sqrt(h1_p) - 30*sqrt(h2_p))/h2_p^3 - 11.1111111111111/h2_p^(5/2);

    B(1,1) = 1.42857142857143/h1_p;
    B(1,2) = 1.42857142857143/h1_p;
    B(2,1) = 0;
    B(2,2) = 0;

    du1 = @(t) F1(t) - F1_p;
    du2 = @(t) FD(t) - FD_p;

    fSL = @(t,dx) A*dx + B*[du1(t); du2(t)];
    dh0 = hl(i, :)' - [h1_p; h2_p];
    [~, dh_temp] = ode45(fSL, [tl(i), tl(i+1)], dh0, odeset('RelTol',1e-6));
    hl(i+1, :) = dh_temp(end, :) + [h1_p, h2_p];
end


%% ---------------------------------------------------------------
% Rysowanie wyników
figure(1);
subplot(2,1,1);
hold on;
plot(tn, hn(:,1), 'b', 'LineWidth', 1.5); 
plot(tl, hl(:,1), 'r--', 'LineWidth', 1.5);
legend('Model nieliniowy', 'Model zlinearizowany', 'Location', 'northeast');
grid on
grid minor
subplot(2,1,2);
hold on;
plot(tn, hn(:,2), 'b', 'LineWidth', 1.5); 
plot(tl, hl(:,2), 'r--', 'LineWidth', 1.5);
legend('Model nieliniowy', 'Model SL', 'Location', 'northeast');
grid on;
grid minor
xlabel('Czas [s]');
ylabel('Poziom cieczy [cm]');


