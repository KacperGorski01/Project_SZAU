%% Sztuczna inteligencja w automatyce - Zadanie 2
% I. Symulacja procesu
clear all;
grubosc = 1.25; set(groot,'DefaultLineLineWidth',grubosc); set(groot,'DefaultStairLineWidth',grubosc);
colors = lines; set(groot, 'defaultAxesColorOrder', colors);

%% 1. Wyznaczanie charakterystyki statycznej y(u)
% Zakres u od u_min do u_max 
u_min = -1; 
u_max = 1;
u_stat_steps = -1:0.05:1; % Próbkowanie zakresu sterowania
y_stat = zeros(size(u_stat_steps));

for i = 1:length(u_stat_steps)
    u_val = u_stat_steps(i);
    % Symulacja d³ugiego kroku, aby osi±gn±æ stan ustalony
    k_sim = 100; 
    u_sim = u_val * ones(1, k_sim);
    y_sim = zeros(1, k_sim);
    x_sim = zeros(1, k_sim);
    
    for k = 6:k_sim % Start od 6 ze wzglêdu na opó¼nienia u(k-5)
        [y_sim(k), x_sim(k)] = proces2_symulator(u_sim(k-4), u_sim(k-5), x_sim(k-1), x_sim(k-2));
    end
    y_stat(i) = y_sim(end); % Pobranie warto¶ci ustalonej [cite: 5]
end

figure;
plot(u_stat_steps, y_stat, 'o-');
title('Charakterystyka statyczna procesu y(u)');
xlabel('u'); ylabel('y');
grid on;

%% 2. Generowanie zbiorów danych (ucz±cy i weryfikuj±cy)
% Parametry generowania danych [cite: 7, 8]
N = 4000;              % Liczba próbek na zbiór
step_period = 50;      % Zmiana sygna³u co 50 kroków

% Inicjalizacja struktur danych
sets = {'ucz±cy', 'weryfikuj±cy'};
data = struct();

for s = 1:2
    u = zeros(1, N);
    y = zeros(1, N);
    x = zeros(1, N);
    
    % Generowanie losowych zmian skokowych w zakresie umin...umax [cite: 6]
    for k = 1:step_period:N
        u_rand = u_min + (u_max - u_min) * rand();
        end_idx = min(k + step_period - 1, N);
        u(k:end_idx) = u_rand;
    end
    
    % Symulacja procesu [cite: 6]
    for k = 6:N
        [y(k), x(k)] = proces2_symulator(u(k-4), u(k-5), x(k-1), x(k-2));
    end
    
    % Zapis do struktury
    data(s).u = u;
    data(s).y = y;
    
    % Wy¶wietlanie danych [cite: 8]
    figure;
    subplot(2,1,1);
    stairs(u);
    title(['Zbiór ', sets{s}, ' - sygna³ steruj±cy u(k)']);
    ylabel('u'); grid on;
    subplot(2,1,2);
    plot(y);
    title(['Zbiór ', sets{s}, ' - wyj¶cie procesu y(k)']);
    ylabel('y'); xlabel('k'); grid on;
end

% Eksport danych do dalszych zadañ (opcjonalnie)
u_train = data(1).u; y_train = data(1).y;
u_verify = data(2).u; y_verify = data(2).y;

% Zapis danych do pliku, aby drugi skrypt móg³ je wczytaæ
save('dane_procesu.mat', 'u_train', 'y_train', 'u_verify', 'y_verify');
disp('Dane zosta³y zapisane do pliku dane_procesu.mat');