%% Porównanie modelu nieliniowego i zlinearizowanego
clear; clc;

Tend = 20000;

% Nieliniowy model dynamiczny

function f = fNonlinear(t,h)
    tau = 125;
    F1 = 200 + 150*(t-tau >= 1000) - 300*(t-tau >= 10000);
    FD = 100 + 100*(t >= 6000) - 200*(t >= 15000);
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
    A = [-0.00740320105820106, 0; 0.0000653086419753086, -0.000111111111111111];
    B = [0.00839682539682540, 0.00839682539682540; 0, 0];
    u1 = 150*(t-tau >= 1000) - 300*(t-tau >= 10000);
    u2 = 100*(t >= 6000) - 200*(t >= 15000);
    f = A*x + B*[u1; u2];
end

x0 = [0; 0];
[t2,x] = ode45(@fLinear, [0, Tend], x0, odeset('RelTol',1e-6));
h1 = x(:,1) + 170.1323;
h2 = x(:,2) + 100;


figure(1);
sgtitle({
    '$F_{in} = 200 + 150 \cdot (t \ge 1000) - 300 \cdot (t \ge 10000)$', ...
    '$F_{D} = 100 + 100 \cdot (t \ge 6000) - 200 \cdot (t \ge 15000)$'
    }, 'Interpreter','latex','FontSize',14);
subplot(2,1,1);
hold on;
plot(t, h(:,1), 'b', 'LineWidth', 1.5); 
plot(t2, h1, 'r--', 'LineWidth', 1.5);
legend('Model nieliniowy', 'Model zlinearizowany', 'Location', 'northeast');
grid on
grid minor
subplot(2,1,2);
hold on;
plot(t, h(:,2), 'b', 'LineWidth', 1.5); 
plot(t2, h2, 'r--', 'LineWidth', 1.5);
legend('Model nieliniowy', 'Model zlinearizowany', 'Location', 'northeast');
grid on;
grid minor
xlabel('Czas [s]');
ylabel('Poziom cieczy [cm]');

% Łatwo zauważyć, że h1 reaguje znacznie szybciej niż h2.
% Bieguny: −0.0074032, −0.00011111


%% Konwencjonalny regulator DMC
clear; clc;

% Punkt linearyzacji:
Fin_point = 200; 
Fd_point = 100;
h1_point = ( 1/23 * (Fin_point + Fd_point) )^2;
h2_point = ( 1/30 * (Fin_point + Fd_point) )^2;

% ---------------------------------------------------------------
% Identyfikacja skoku jednostkowego dla sterowania 

Ts = 150; % okres próbkowania dla regulatora DMC
Tend = 50000;
T = Ts : Ts : Tend;

tau = 125;
A = [-0.00740320105820106, 0; 0.0000653086419753086, -0.000111111111111111];
B = [0.00839682539682540, 0.00839682539682540; 0, 0];

f1 = @(t_,x_) A*x_ + B*[(t_ >= 0); 0];
[t, x] = ode45(f1, [0, Tend], [0; 0], odeset('RelTol',1e-6));
y1 = x(:,2);
s1 = interp1(t, y1, T);

f2 = @(t_,x_) A*x_ + B*[0; (t_ >= 0)];
[t, x] = ode45(f2, [0, Tend], [0; 0], odeset('RelTol',1e-6));
y2 = x(:,2);
s2 = interp1(t, y2, T);

figure(1)
hold on
plot(t, y1, '-', T, s1, 'o')
plot(t, y2, '-', T, s2, 'o')
title('Odpowiedzi skokowe modelu zlinearizowanego')
xlabel('Czas [s]')
grid on
grid minor

% Łatwo zauważyć, że odpowiedzi skokowe dla obu wejść są zawsze identyczne - i to nie zależnie od punktu w jakim linearyzujemy. 
% Więc możemy przyjąć s = s1 = s2. (Macierze dynamiczne będą takie same)
s = s1;

clear t x y1 y2 s1 s2 T

% ---------------------------------------------------------------
% Wyznaczenie macierzy DMC

D = length(s);      % horyzont dynamiki
N = round(D/100);   % horyzont predykcji
Nu = 10;            % horyzont sterowania
lambda = 4;         % współczynnik kary

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
K1 = K(1,:);  
ke = sum(K1);  
kp = K1 * Mp; 

% ---------------------------------------------------------------
% Symulacja działania regulatora DMC na modelu nieliniowym
% Matlab numeruje indeksy od 1 !!!
% Czyli dla indeksu 1 mamy chwilę '0*Ts sek', dla indeksu 2 mamy chwilę '1*Ts sek', itd.

time = 0 : Ts : 1800000;             % czas symulacji

h = [170, 100];  % początkowa wartość stanów

y = zeros(length(time), 1);     % wektor wyjść
y(1) = h(1,2);                  % początkowa wartość wyjścia   

y_zad = 100 + 20 * (time >= 400000) - 20 * (time >= 1000000);   % wartość zadana
Fd = 100 + 100 * (time >= 700000) - 100 * (time >= 1300000);    % zakłócenia

Fin = zeros(length(time), 1);     % wektor sterowań

dUp1 = zeros(D-1, 1);            % wektor przyrostów sterowania z poprzednich kroków (początkowo zerowy)
dUp2 = zeros(D-1, 1);            % wektor przyrostów zakłócenia z poprzednich kroków (początkowo zerowy)


for k = 1 : length(time) - 1
    % Aktualizacja wektora przyrostów zakłócenia z poprzednich kroków
    if k > 1
        dUp2 = [Fd(k) - Fd(k-1); dUp2(1 : end - 1)]; 
    end

    % Obliczanie przyrostu sterowania
    du = ke * ( y_zad(k) - y(k) ) - dot(kp, dUp1 + dUp2);
    % Aktualizacja wektora przyrostów sterowania z poprzednich kroków
    dUp1 = [du; dUp1(1 : end - 1)];
    
    % Aktualizacja sterowania
    if k > 1
        Fin(k) = Fin(k - 1) + du;
    else
        Fin(k) = du + Fin_point; % dla k=1 mamy przyrost sterowania + steroanie z punktu równowagi
    end
    
    % Przy nowym sterowaniu obliczamy wyjście w następnej chwili próbkowania.
    f2 = @(t_,h_) [
        ( Fin(k) + Fd(k) - 23 * sqrt(h_(1))) / (0.7 * (h_(1))) ;
        ( 23 * sqrt(h_(1)) - 30 * sqrt(h_(2))) / (1.35 * (h_(2))^2)
    ];
    [~, h_temp] = ode45(f2, [time(k) time(k+1)], h, odeset('RelTol',1e-3));
    h = h_temp(end,:);
    y(k+1) = h(2);
end

% rysujemy wykresy
% (poprawiamy jeszcze ostatnią wartość sterowania, aby wykres był ładniejszy (zabieg kosmetyczny))
Fin(end) = Fin(end-1);
figure(2)
subplot(3,1,1)
hold on
plot(time, y, 'LineWidth', 1.5)
plot(time, y_zad, 'LineWidth', 1.5)
ylim([min(y)-5, max(y)+5])
grid on
grid minor
title('Wyjście y i wartość zadana [cm]')
xlabel('Czas [s]')
subplot(3,1,2)
plot(time, Fin, 'LineWidth', 1.5)
ylim([min(Fin)-5, max(Fin)+5])
grid on
grid minor
title('Sterowanie F_{in} [cm^3/s]')
xlabel('Czas [s]')
subplot(3,1,3)
plot(time, Fd, 'LineWidth', 1.5)
ylim([min(Fd)-5, max(Fd)+5])
title('Zakłócenie F_{D} [cm^3/s]')
xlabel('Czas [s]')
grid on
grid minor


%% 
clear; clc;





