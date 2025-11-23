%% Zwykły liniowy regulator DMC
clear; clc;

Ts = 125/2; % okres próbkowania dla regulatora DMC
% dla wejścia bez opóźnienia działało dobrze Ts = 100 sek, dlatego uwzględniając opóźnienie 125s przyjeliśmy 125/2 sek
% tak, że opóźnienie wynosi dokładnie 2 okresy próbkowania (łatwiej wtedy implementować opóźnienie w symulacji)

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

Tend = 50000;
TT = Ts : Ts : Tend;

tau = 125;
f1 = @(t_,x_) A*x_ + B*[(t_ - tau >= 0); 0];
[t1, x1] = ode45(f1, [0, Tend], [0; 0], odeset('RelTol',1e-6));
y1 = x1(:,2);
s1 = interp1(t1, y1, TT);

f2 = @(t_,x_) A*x_ + B*[0; (t_ - tau >= 0)];
[t2, x2] = ode45(f2, [0, Tend], [0; 0], odeset('RelTol',1e-6));
y2 = x2(:,2);
s2 = interp1(t2, y2, TT);

figure(1)
hold on
plot(t1, y1, '-', TT, s1, 'o')
plot(t2, y2, '-', TT, s2, 'o')
title('Odpowiedzi skokowe modelu zlinearizowanego')
xlabel('Czas [s]')
grid on
grid minor

% Łatwo zauważyć, że odpowiedzi skokowe dla obu wejść są zawsze identyczne - i to nie zależnie od punktu w jakim linearyzujemy.
% Wynika to z równości odpowiednich pochodnych funkcji nieliniowej opisującej układ. 
% Więc możemy przyjąć s = s1 = s2. (Macierze dynamiczne będą takie same)
s = s1;

clear t1 t2 x1 x2 y1 y2 s1 s2 TT tau f1 f2


%% ---------------------------------------------------------------
% Wyznaczenie parametrów regulatora

D = length(s);      % horyzont dynamiki
N = round(D/100);   % horyzont predykcji
Nu = 10;            % horyzont sterowania
lambda = 10;        % współczynnik kary

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

Tend = 12e5;

time = 0 : Ts : Tend;             % czas symulacji

h = [170, 100];  % początkowa wartość stanów

y = zeros(length(time), 1);     % wektor wyjść (h2)
y(1) = h(1,2);                  % początkowa wartość wyjścia   

y_zad = 100 + 150*(time>0.1e5) - 200*(time>6e5);   % wartość zadana
Fd = 100 + 50*(time>2.5e5) - 100*(time>8.5e5);    % zakłócenia

Fin = zeros(length(time), 1);     % wektor sterowań

dUp1 = zeros(D-1, 1);            % wektor przyrostów sterowania z poprzednich kroków (początkowo zerowy)
dUp2 = zeros(D-1, 1);            % wektor przyrostów zakłócenia z poprzednich kroków (początkowo zerowy)

for k = 1 : length(time) - 1
    fprintf('[%f %%]\n', 100 * k / (length(time)-1));

    % Obliczanie przyrostu sterowania
    du1 = ke * ( y_zad(k) - y(k) ) - dot(kp, dUp1 + dUp2);
    du1 = min(max(du1, -50), 50);  % ograniczenie przyrostu sterowania

    % Aktualizacja wektora przyrostów sterowania z poprzednich kroków
    dUp1 = [du1; dUp1(1 : end - 1)];

    % Aktualizacja wektora przyrostów zakłócenia z poprzednich kroków
    if k > 1
        dUp2 = [Fd(k) - Fd(k-1); dUp2(1 : end - 1)]; 
    end
    
    % Aktualizacja sterowania
    if k > 1
        Fin(k) = Fin(k - 1) + du1;
    else
        Fin(k) = du1 + Fin_point; % dla k=1 mamy przyrost sterowania + steroanie z punktu równowagi
    end

    Fin(k) = min(max(Fin(k), 0), 700); % ograniczenie sterowania
    
    % Przy nowym sterowaniu obliczamy wyjście w następnej chwili próbkowania.
    % Uwzględniamy opóźnienie sterowania równe 125s (czyli 2 okresy próbkowania Ts)
    if k > 2
        Fin_ = Fin(k-2);
    else
        Fin_ = Fin_point;
    end
    f2 = @(t_,h_) [ ...
        ( Fin_ + Fd(k) - 23 * sqrt(max(h_(1), 1e-6))) / (0.7 * max(h_(1), 1e-6)) ; ...
        ( 23 * sqrt(max(h_(1), 1e-6)) - 30 * sqrt(max(h_(2), 1e-6))) / (1.35 * (max(h_(2), 1e-6))^2) ...
    ];
    [~, h_temp] = ode45(f2, [time(k) time(k+1)], h, odeset('RelTol',1e-3));
    h = h_temp(end,:);
    h(1) = max(h(1), 1e-6);
    h(2) = max(h(2), 1e-6);
    y(k+1) = h(2);
end


%% ---------------------------------------------------------------
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



