%% Rozmyty DMC
clear; clc;

Ts = 125/2; % okres próbkowania dla regulatora DMC (jest wspólny dla wszystkich modeli)
% dla wejścia bez opóźnienia działało dobrze Ts = 100 sek, dlatego uwzględniając opóźnienie 125s przyjeliśmy 125/2 sek
% tak, że opóźnienie wynosi dokładnie 2 okresy próbkowania (łatwiej wtedy implementować opóźnienie w symulacji)


%% ---------------------------------------------------------------
% Identyfikacja skoku jednostkowego dla sterowania
% i przygotowanie parametrów regulatorów DMC.

% Zbiory rozmyte dla wyjścia h2
h2Range = [100, 500];

NumOfFuzzySets = 5;
h2p = linspace(h2Range(1), h2Range(2), NumOfFuzzySets);
Dh2p = h2p(2) - h2p(1);

mf = cell(1, NumOfFuzzySets);
sigma = 0.9*Dh2p;
for i = 1:NumOfFuzzySets
    c = h2p(i);
    mf{i} = @(z) gaussmf(z, [sigma c]);
end

% Lokalne modele zlinearizowane

Tend = [5e4, 2e5, 5e5, 1e6, 2e6];
T = cell(1, NumOfFuzzySets);
for i = 1:NumOfFuzzySets
    T{i} = Ts : Ts : Tend(i);
end

% W celu weryfikacji wyświetlaliśmy pobrane odpowiedzi skokowe (teraz zakomentowane)
% figure(1);
% hold on

s = cell(1, NumOfFuzzySets);
Mp = cell(1, NumOfFuzzySets);
M = cell(1, NumOfFuzzySets);
K = cell(1, NumOfFuzzySets);

ke = cell(1, NumOfFuzzySets);
kp = cell(1, NumOfFuzzySets);

D = zeros(1, NumOfFuzzySets);
N = zeros(1, NumOfFuzzySets);

for i = 1:NumOfFuzzySets
    h2p_ = h2p(i);
    Fdp_ = 100;
    F1p_ = sqrt(h2p_)*30 - Fdp_;
    h1p_ = ( 1/23 * (F1p_ + Fdp_) )^2;

    A = zeros(2,2);
    A(1,1) = -1.42857142857143*(F1p_ + Fdp_ - 23*sqrt(h1p_))/h1p_^2 - 16.4285714285714/h1p_^(3/2);
    A(1,2) = 0;
    A(2,1) = 8.51851851851852/(sqrt(h1p_)*h2p_^2);
    A(2,2) = -1.48148148148148*(23*sqrt(h1p_) - 30*sqrt(h2p_))/h2p_^3 - 11.1111111111111/h2p_^(5/2);

    B = zeros(2,2);
    B(1,1) = 1.42857142857143/h1p_;
    B(1,2) = 1.42857142857143/h1p_;
    B(2,1) = 0;
    B(2,2) = 0;

    fLin = @(t_,x_) A*x_ + B*[(t_ - 125 >= 0); 0];
    [t, x] = ode45(fLin, [0, Tend(i)], [0; 0], odeset('RelTol',1e-6));
    y = x(:,2);
    s{i} = interp1(t, y, T{i});

    % plot(t, y, '-', T{i}, s{i}, 'o')

    D(i) = length(s{i});      % horyzont dynamiki
    N(i) = round(D(i)/100);   % horyzont predykcji
    Nu = 10;                  % horyzont sterowania
    lambda = 10;              % współczynnik kary
    
    Mp{i} = zeros(N(i), D(i)-1);
    for r = 1 : N(i)
    for j = 1 : D(i)-1
        if r + j <= D(i)
            Mp{i}(r,j) = s{i}(r + j) - s{i}(j);
        else
            Mp{i}(r,j) = s{i}(D(i)) - s{i}(j);
        end
    end
    end
    
    M{i} = zeros(N(i), Nu);
    for r = 1 : N(i)
    for j = 1 : Nu
        if r >= j
            M{i}(r,j) = s{i}(r - j + 1);
        else
            M{i}(r,j) = 0;
        end
    end
    end
    
    K{i} = (M{i}' * M{i} + lambda * eye(Nu)) \ M{i}';

    ke{i} = sum(K{i}(1,:));  
    kp{i} = K{i}(1,:) * Mp{i}; 
end
% grid on
% grid minor

clear M Mp K t x y s


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

dUp1 = zeros(max(D)-1, 1);            % wektor przyrostów sterowania z poprzednich kroków (początkowo zerowy)
dUp2 = zeros(max(D)-1, 1);            % wektor przyrostów zakłócenia z poprzednich kroków (początkowo zerowy)

for k = 1 : length(time) - 1
    fprintf('[%f %%]\n', 100 * k / (length(time)-1));

    % Obliczanie przyrostu sterowania
    du1__ = zeros(1, NumOfFuzzySets);
    mf__ = zeros(1, NumOfFuzzySets);
    for j = 1 : NumOfFuzzySets
        mf__(j) = mf{j}(y(k));
        du1__ (j) = ke{j} * ( y_zad(k) - y(k) ) - dot(kp{j}, dUp1(1 : D(j)-1) + dUp2(1 : D(j)-1));
    end
    du1 = dot(mf__, du1__) / sum(mf__);
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










