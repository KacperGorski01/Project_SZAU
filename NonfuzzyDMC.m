%% Konwencjonalny regulator DMC
clear; clc;


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
% Identyfikacja skoku jednostkowego dla sterowania 

Ts = 150; % okres próbkowania dla regulatora DMC
Tend = 50000;
T = Ts : Ts : Tend;

tau = 125;
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


%% ---------------------------------------------------------------
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

clear Mp M K K1


%% ---------------------------------------------------------------
% Symulacja działania regulatora DMC na modelu nieliniowym
% Matlab numeruje indeksy od 1 !!!
% Czyli dla indeksu 1 mamy chwilę '0*Ts sek', dla indeksu 2 mamy chwilę '1*Ts sek', itd.

time = 0 : Ts : 1800000;             % czas symulacji

h = [170, 100];  % początkowa wartość stanów

y = zeros(length(time), 1);     % wektor wyjść
y(1) = h(1,2);                  % początkowa wartość wyjścia   

y_zad = 100 + 20 * (time >= 400000) - 20 * (time >= 1000000);   % wartość zadana
Fd = 100 + 10 * (time >= 700000) - 10 * (time >= 1300000);    % zakłócenia

Fin = zeros(length(time), 1);     % wektor sterowań

dUp1 = zeros(D-1, 1);            % wektor przyrostów sterowania z poprzednich kroków (początkowo zerowy)
dUp2 = zeros(D-1, 1);            % wektor przyrostów zakłócenia z poprzednich kroków (początkowo zerowy)

for k = 1 : length(time) - 1
    % Aktualizacja wektora przyrostów zakłócenia z poprzednich kroków
    if k > 1
        dUp2 = [Fd(k) - Fd(k-1); dUp2(1 : end - 1)]; 
    end

    % Obliczanie przyrostu sterowania
    du1 = ke * ( y_zad(k) - y(k) ) - dot(kp, dUp1 + dUp2);
    % Aktualizacja wektora przyrostów sterowania z poprzednich kroków
    dUp1 = [du1; dUp1(1 : end - 1)];
    
    % Aktualizacja sterowania
    if k > 1
        Fin(k) = Fin(k - 1) + du1;
    else
        Fin(k) = du1 + Fin_point; % dla k=1 mamy przyrost sterowania + steroanie z punktu równowagi
    end
    
    % Przy nowym sterowaniu obliczamy wyjście w następnej chwili próbkowania.
    f2 = @(t_,h_) [ ...
        ( Fin(k) + Fd(k) - 23 * sqrt(h_(1))) / (0.7 * (h_(1))) ; ...
        ( 23 * sqrt(h_(1)) - 30 * sqrt(h_(2))) / (1.35 * (h_(2))^2) ...
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



