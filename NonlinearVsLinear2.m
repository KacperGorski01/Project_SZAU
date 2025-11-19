%% Porównanie modelu nieliniowego i zlinearizowanego
clear; clc;

Tend = 2e6;

% Uwaga:
% Można zauważyć, że jakość dopasowania modelu zlinearyzowanego do
% nieliniowego nie zależy osobno od F1 i Fd, ale od ich sumy F1+Fd.
% Jest to zresztą zgodne z fizyką zagadnienia.
% Przybliżenie działa w miarę sensownie o ile wychylenie od 
% sumy (z punktu linearyzacji-równowagi) nie przekracza +/-20%.


%% ---------------------------------------------------------------
% Linearyzacja modelu nieliniowego

% Punkt linearyzacji (punkt równowagi):
Fin_point = 200; % (F1_point)
Fd_point = 100;
h1_point = ( 1/23 * (Fin_point + Fd_point) )^2;
h2_point = ( 1/30 * (Fin_point + Fd_point) )^2;

A = zeros(2,2);
A(1,1) = -1.42857142857143*(Fin_point + Fd_point - 23*sqrt(h1_point))/h1_point^2 - 16.4285714285714/h1_point^(3/2);
A(1,2) = 0;
A(2,1) = 8.51851851851852/(sqrt(h1_point)*h2_point^2);
A(2,2) = -1.48148148148148*(23*sqrt(h1_point) - 30*sqrt(h2_point))/h2_point^3 - 11.1111111111111/h2_point^(5/2);

B = zeros(2,2);
B(1,1) = 1.42857142857143/h1_point;
B(1,2) = 1.42857142857143/h1_point;
B(2,1) = 0;
B(2,2) = 0;


%% ---------------------------------------------------------------
% Symulacja modeli nieliniowego i zlinearizowanego

% Wychylenia od punktu równowagi:
dF1 = @(t) 100*(t>= 2e5) - 100*(t >= 1e6);
dFD = @(t) 50*(t >= 5e5) - 50*(t >= 1.5e6);

% Nieliniowy model dynamiczny
h0 = [ h1_point ; h2_point ];
F1 = @(t) Fin_point + dF1(t);
FD = @(t) Fd_point + dFD(t);

fNonlinear = @(t,h) [ ...
    ( F1(t) + FD(t) - 23 * sqrt(h(1))) / (0.7 * h(1)) ; ...
    ( 23 * sqrt(h(1)) - 30 * sqrt(h(2))) / (1.35 * h(2)^2) ...
];

[tn,hn] = ode45(fNonlinear, [0, Tend], h0, odeset('RelTol',1e-6));

% Zlinearizowany model dynamiczny

% h1 = h1_point + dh1
% h2 = h2_point + dh2
% F1 = F1_point + dF1
% FD = FD_point + dFD

dh0 = [0; 0];
fLinear = @(t,dh) A*dh + B*[dF1(t); dFD(t)];
[tl, dhl] = ode45(fLinear, [0, Tend], dh0, odeset('RelTol',1e-6));
hl = dhl + [h1_point, h2_point];

% Zlinearyzowany model dynamiczny - rozwiązany metodą Eulera

Ts = Tend/5e4;
te = 0:Ts:Tend;

he = zeros(2, length(te));
he(:,1) = h0;
for k = 1:length(te)-1
    dh = A*(he(:,k) - h0) + B*([F1(te(k)); FD(te(k))] - [Fin_point; Fd_point]);
    he(:,k+1) = he(:,k) + dh*Ts;
end


% Rysowanie wyników
figure(1);
subplot(2,1,1);
hold on;
plot(tn, hn(:,1), 'b', 'LineWidth', 1.5); 
plot(tl, hl(:,1), 'r--', 'LineWidth', 1.5);
plot(te, he(1,:), 'g:', 'LineWidth', 1.5);
legend('Model nieliniowy', 'Model zlinearizowany', 'Model zlinearizowany (Euler)', 'Location', 'northeast');
grid on
grid minor
subplot(2,1,2);
hold on;
plot(tn, hn(:,2), 'b', 'LineWidth', 1.5); 
plot(tl, hl(:,2), 'r--', 'LineWidth', 1.5);
plot(te, he(2,:), 'g:', 'LineWidth', 1.5);
legend('Model nieliniowy', 'Model zlinearizowany', 'Model zlinearizowany (Euler)', 'Location', 'northeast');
grid on;
grid minor


