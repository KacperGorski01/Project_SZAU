% --- Parametry DMC ---
Nu = 10;            % Horyzont sterowania
N = 200;            % Horyzont predykcji
Du = 1600;          
lambda = 2;         % Lambda

ks = 1500;          % Początek sterowania
kp = 1400;          % Punkt przejściowy
kk = 8000;          % Liczba kroków symulacji
T = 1;              % Krok dyskretyzacji

% --- Punkt linearyzacji ---
C1 = 0.75;          % Stała dla zbiornika 1
C2 = 0.55;          % Stała dla zbiornika 2
alpha1 = 20;        % Współczynnik przepływu z 1 do 2
alpha2 = 20;        % Współczynnik odpływu

% --- Symulacja odpowiedzi skokowej ---
F1_0 = 52;          % Skok wejściowy
FD_0 = 11;          % Zakłócenie
h1_0 = 9.9225;      % Początkowy poziom h1
h2_0 = 5.9225;      % Początkowy poziom h2
v1_0 = h1_0^2 * C1; % Objętość h1
v2_0 = h2_0^3 * C2; % Objętość h2

h1 = h1_0 * ones(1, kk);
h2 = h2_0 * ones(1, kk);
v1 = h1_0^2 * C1 * ones(1, kk);
v2 = h2_0^3 * C2 * ones(1, kk);
h2_step = h1_0;
F1_step = alpha2 * sqrt(h2_step) - FD_0;
F1in = F1_0 * ones(1, kk);

step_time = 0:1:kk-1; % Czas symulacji odpowiedzi skokowej
h2_step_values = zeros(size(step_time));
dt_step = step_time(2) - step_time(1);

% Pętla symulacji odpowiedzi skokowej
for k = 2:kk
    if k/T > 1
        F1in(k) = F1_step;
    end

    v1(k) = v1(k-1) + T * (F1in(k-1) - F1_0 + FD_0 - (alpha1 / (2 * sqrt(h1_0))) * (h1(k-1) - h1_0));
    v2(k) = v2(k-1) + T * ((alpha1 / (2 * sqrt(h1_0))) * (h1(k-1) - h1_0) - (alpha2 / (2 * sqrt(h2_0))) * (h2(k-1) - h2_0));
    h2(k) = h2_0 + 1 / (3 * (v2_0^2 * C2)^(1/3)) * (v2(k) - v2_0);
    h1(k) = h1_0 + 1 / (2 * sqrt(v1_0 * C1)) * (v1(k) - v1_0);
end

% Obliczanie odpowiedzi skokowej systemu
step_response = (h2 - h2_0) / h2_step;

% Regulator DMC
M = zeros(N, Nu);
for i = 1:N
    for j = 1:Nu
        if (i >= j)
            M(i, j) = step_response(i - j + 1);
        end
    end
end

Mp = zeros(N, Du - 1);
for i = 1:N
    for j = 1:Du-1
        if i + j <= Du
            Mp(i, j) = step_response(i + j) - step_response(j);
        else
            Mp(i, j) = step_response(Du) - step_response(j);
        end
    end
end

% Wyznaczanie parametrów prawa sterowania
I = eye(Nu);
K = ((M' * M + lambda * I) ^ -1) * M';
ke = sum(K(1, :));
ku = K(1, :) * Mp;
deltaup = zeros(1, Du-1);

% Symulacja systemu z DMC
h1_dmc = h1_0;
h2_dmc = h2_0;
v1_dmc = v1_0;
v2_dmc = v2_0;
h2_dmc_values = zeros(1, numel(step_time));
delta_u = zeros(1, numel(step_time));
h2_zad = h2_0 * ones(1, kk);
h2_zad(kp:kk) = h1_0;
F1 = F1_0;
FD = FD_0;

for k = 2:numel(step_time)
    y = h2_dmc; % Odchylenie od punktu pracy
    delta_u(k) = ke * (h2_zad(k) - y) - ku * deltaup';
    deltaup = circshift(deltaup, 1);
    deltaup(1) = delta_u(k);

    F1 = F1 + delta_u(k); % Aktualizacja sygnału sterowania
    v1_dmc = v1_dmc + T * (F1 + FD_0 - alpha1 * sqrt(h1_dmc));
    v2_dmc = v2_dmc + T * (alpha1 * sqrt(h1_dmc) - alpha2 * sqrt(h2_dmc));
    h1_dmc = sqrt(v1_dmc / C1);
    h2_dmc = (v2_dmc / C2)^(1/3);
    h2_dmc_values(k) = h2_dmc;
end

% --- Rysowanie wykresów ---
figure('Position', [100, 50, 800, 400]);
plot(step_time, h2_zad, 'b');
hold on;
plot(step_time, h2_dmc_values(1:k), 'r');
xlabel('Czas [s]');
ylabel('h2 [cm]');
title('Regulator DMC analityczny');
legend('Wartość zadana h2', 'Wyjście h2');
grid on;
print(gcf, '1.16.png', '-dpng', '-r300', '-opengl');
