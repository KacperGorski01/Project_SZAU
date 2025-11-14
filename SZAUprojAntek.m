%% Porównanie modelu nieliniowego i zlinearizowanego
clear; clc;

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


%% Konwencjonalny regulator DMC
clear; clc;

% x1 = delta_h1 = h1 - 170.1323
% x2 = delta_h2 = h2 - 100
% u1 = delta_F1 = F1 - 200
% u2 = delta_FD = FD - 100
% y = x2 = h2 - 100

Ts = 100; % okres próbkowania dla regulatora DMC

% ---------------------------------------------------------------
% Identyfikacja skoku jednostkowego dla sterowania 

tau = 125;
A = [-0.0074032, 0; 0.0000653, -0.00011111];
B = [0.00839683; 0];

% Bieguny: −0.0074032, −0.00011111
% wyjście y jest bardzo wolne

f1 = @(t_,x_) A*x_ + B*(t_ >= 0);
[t, x] = ode45(f1, [0, 70000], [0; 0], odeset('RelTol',1e-6));
y = x(:,2);

T = Ts : Ts : 60000;        % wektor czasu próbkowania
s = interp1(t, y, T);      % wektor odpowiedzi skokowej w chwilach próbkowania

figure(1)
hold on
plot(t, y, '-', T, s, 'o')
grid on
grid minor

% ---------------------------------------------------------------
% Wyznaczenie macierzy DMC

D = length(s);      % horyzont dynamiki
N = round(D/100);     % horyzont predykcji
Nu = 10;             % horyzont sterowania
lambda = 10;         % współczynnik kary

Mp = zeros(N, D-1);
for i = 1 : N
for j = 1 : D-1
    if i + j <= D
        Mp(i,j) = s(i + j) - s(j);
    else
        Mp(i,j) = s(D) - s(j);
    end
end
end

M = zeros(N, Nu);
for i = 1 : N
for j = 1 : Nu
    if i >= j
        M(i,j) = s(i - j + 1);
    else
        M(i,j) = 0;
    end
end
end

K = (M' * M + lambda * eye(Nu)) \ M';

ke = sum(K(1,:));
ku = K(1,:) * Mp;

% ---------------------------------------------------------------
% Symulacja działania regulatora DMC na modelu nieliniowym

% Matlab numeruje indeksy od 1 !!!
% Czyli dla indeksu 1 mamy chwilę '0*Ts sek', dla indeksu 2 mamy chwilę '1*Ts sek', itd.

time = 0 : Ts : 2000000;             % wektor czasu symulacji

x = zeros(length(time), 2);     % wektor stanów
x(1,:) = [0; 0];  % początkowa wartość stanów

y = zeros(length(time), 1);     % wektor wyjść
y(1) = x(1,2);                      % początkowa wartość wyjścia   

y_zad = 5 * (time >= 500000) + 5 * (time >= 1000000) - 5 * (time >= 1500000);        % wartość zadana

u = zeros(length(time), 1);     % wektor sterowań

dUp = zeros(D-1, 1);            % wektor przyrostów sterowania z poprzednich kroków (początkowo zerowy)

for k = 1 : length(time) - 1
    % Obliczanie przyrostu sterowania
    du = ke * (y_zad(k) - y(k)) - dot(ku, dUp);
    
    % Aktualizacja wektora przyrostów sterowania z poprzednich kroków
    dUp = [du; dUp(1 : end - 1)];
    
    % Aktualizacja sterowania
    if k > 1
        u(k) = u(k - 1) + du;
    else
        u(k) = du; % dla k=1
    end
    
    % Przy nowym sterowaniu obliczamy wyjście w następnej chwili próbkowania.
    f2 = @(t_,x_) A*x_ + B*u(k);
    [~, x_temp] = ode45(f2, [time(k) time(k+1)], x(k,:), odeset('RelTol',1e-3));
    x(k+1,:) = x_temp(end,:);
    y(k+1) = x(k+1,2);
end

figure(2)
hold on
plot(time, y)
plot(time, y_zad)
grid on
grid minor

figure(3)
hold on
stairs(time, u)
grid on
grid minor



