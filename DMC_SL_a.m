%% DMC z sukcesywną linearyzacją - analityczny

% Uwaga: Z racji na dużo obliczeń w każdej iteracji symulacji skrypt może działać dość długo (nawet 30min)

clear; clc;

Tend = 12e5;
Ts = 125/2; % okres próbkowania dla regulatora DMC
% dla wejścia bez opóźnienia działało dobrze Ts = 100 sek, dlatego uwzględniając opóźnienie 125s przyjeliśmy 125/2 sek
% tak, że opóźnienie wynosi dokładnie 2 okresy próbkowania (łatwiej wtedy implementować opóźnienie w symulacji)

time = 0 : Ts : Tend;             % czas symulacji

h = [170, 100];  % początkowa wartość stanów

y = zeros(length(time), 1);     % wektor wyjść (h2)
y(1) = h(1,2);                  % początkowa wartość wyjścia   

y_zad = 100 + 150*(time>0.1e5) - 200*(time>6e5);   % wartość zadana
Fd = 100 + 50*(time>2.5e5) - 100*(time>8.5e5);    % zakłócenia

Fin = zeros(length(time), 1);     % wektor sterowań

D = 15e3;           % horyzont dynamiki
N = 500;            % horyzont predykcji
Nu = 10;            % horyzont sterowania
lambda = 10;        % współczynnik kary

TT = Ts:Ts:D*Ts;

s = zeros(1, D);
Mp = zeros(N, D-1);
M = zeros(N, Nu);
K = zeros(Nu, N);
ke = 0;
kp = zeros(1, D-1);

dUp1 = zeros(D-1, 1);  % wektor przyrostów sterowania z poprzednich kroków (początkowo zerowy)
dUp2 = zeros(D-1, 1);  % wektor przyrostów zakłócenia z poprzednich kroków (początkowo zerowy)

A = zeros(2,2);
B = zeros(2,2);
for k = 1 : length(time) - 1
    fprintf('[%f %%]\n', 100 * k / (length(time)-1));

    % Linearyzacja i wyznaczenie współczynników DMC
    h2_p = y(k);
    Fd_p = Fd(k);
    h1_p = (30/23)^2 * h2_p;
    F1_p = 23 * sqrt(h1_p) - Fd_p;

    A(1,1) = -1.42857142857143*(F1_p + Fd_p - 23*sqrt(h1_p))/h1_p^2 - 16.4285714285714/h1_p^(3/2);
    A(1,2) = 0;
    A(2,1) = 8.51851851851852/(sqrt(h1_p)*h2_p^2);
    A(2,2) = -1.48148148148148*(23*sqrt(h1_p) - 30*sqrt(h2_p))/h2_p^3 - 11.1111111111111/h2_p^(5/2);

    B(1,1) = 1.42857142857143/h1_p;
    B(1,2) = 1.42857142857143/h1_p;
    B(2,1) = 0;
    B(2,2) = 0;

    [tTmp, xTmp] = ode45(@(t_,x_) A*x_ + B*[(t_ - 125 >= 0);0], [0;D*Ts], [0;0]);
    s = interp1(tTmp, xTmp(:,2), TT);

    for i = 1 : N
    for j = 1 : D-1
        if i + j <= D
            Mp(i,j) = s(i + j) - s(j);
        else
            Mp(i,j) = s(D) - s(j);
        end
    end
    end

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
    kp = K(1,:) * Mp;

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
        Fin(k) = du1 + 200;
    end

    Fin(k) = min(max(Fin(k), 0), 700); % ograniczenie sterowania
    
    % Przy nowym sterowaniu obliczamy wyjście w następnej chwili próbkowania.
    % Uwzględniamy opóźnienie sterowania równe 125s (czyli 2 okresy próbkowania Ts)
    if k > 2
        Fin_ = Fin(k-2);
    else
        Fin_ = 200;
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
figure(4)
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


