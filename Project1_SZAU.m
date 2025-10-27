% Przedmiot: Sztuczna Inteligencja w Automatyce
% Model: Przepływ cieczy między dwoma zbiornikami
% Autorzy: Kacper Górski 314061, Antoni Marczuk

clearvars;
close all;
clc;

folderPath = fullfile(pwd, 'Wykresy');
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end


% parametry modelu
C1 = 0.35;
C2 = 0.45;
alpha1 = 23;
alpha2 = 30;
tau = 125;

F1_0 = 200;
FD_0 = 100;
F2_0 = F1_0 + FD_0;

h1_0 = (F2_0/alpha1)^2;
h2_0 = (F2_0/alpha2)^2;

fprintf('Punkt pracy:\n h1_0 = %.2f cm, h2_0 = %.2f cm\n', h1_0, h2_0);

% parametry symulacji
T_end = 6500;
dt = 1;
t = 0:dt:T_end;

Fin = F1_0 * ones(size(t));
FD = FD_0 * ones(size(t));

Fin(t > 500  & t <= 2000) = 1.1 * F1_0;
Fin(t > 2000 & t <= 3500) = 0.9 * F1_0;
Fin(t > 3500 & t <= 5000) = 1.2 * F1_0;
Fin(t > 5000) = 0.8 * F1_0;

FD(t > 1000 & t <= 3500) = 120;
FD(t > 3500) = 80;

% opóźnienie
N_delay = round(tau / dt);
Fin_delay = [F1_0 * ones(1, N_delay), Fin(1:end - N_delay)];

model_nieliniowy = @(t, x, Fin_t, FD_t) [
    (Fin_t + FD_t - alpha1 * sqrt(max(x(1), 1e-6))) / (2 * C1 * max(x(1), 1e-6));
    (alpha1 * sqrt(max(x(1), 1e-6)) - alpha2 * sqrt(max(x(2), 1e-6))) / (3 * C2 * max(x(2)^2, 1e-6))
];

% symulacja sterowania
x = [h1_0; h2_0];
h1 = zeros(size(t));
h2 = zeros(size(t));

for k = 1:length(t)
    h1(k) = x(1);
    h2(k) = x(2);
    dx = model_nieliniowy(t(k), x, Fin_delay(k), FD_0);
    x = x + dx * dt;
end

figure('Name', 'Symulacja obiektu nieliniowego - zmiana sterowania', 'Color', 'w');
subplot(3,1,1);
plot(t, Fin, 'k--', 'LineWidth', 1.2); hold on;
plot(t, Fin_delay, 'r', 'LineWidth', 1.5);
xlabel('Czas [s]'); ylabel('F_{in} [cm^3/s]');
title('Sterowanie F_{in}');
legend('Wejście zadane', 'Wejście opóźnione');
grid on;

subplot(3,1,2);
plot(t, h1, 'b', 'LineWidth', 1.5);
xlabel('Czas [s]'); ylabel('Poziom h_1 [cm]');
title('Poziom cieczy w pierwszym zbiorniku');
grid on;

subplot(3,1,3);
plot(t, h2, 'r', 'LineWidth', 1.5);
xlabel('Czas [s]'); ylabel('Poziom h_2 [cm]');
title('Poziom cieczy w drugim zbiorniku');
grid on;
print(fullfile(folderPath, 'Obiekt nieliniowy - zmiana sterowania.png'), '-dpng', '-r400');

% symulacja zakłócenia
x = [h1_0; h2_0];
h1 = zeros(size(t));
h2 = zeros(size(t));

for k = 1:length(t)
    h1(k) = x(1);
    h2(k) = x(2);
    dx = model_nieliniowy(t(k), x, F1_0, FD(k));
    x = x + dx * dt;
end

figure('Name', 'Symulacja obiektu nieliniowego - zmiana zakłócenia', 'Color', 'w');
subplot(3,1,1);
plot(t, FD, 'LineWidth', 1.5);
xlabel('Czas [s]'); ylabel('F_D [cm^3/s]');
title('Zakłócenie F_D(t)');
grid on;

subplot(3,1,2);
plot(t, h1, 'b', 'LineWidth', 1.5);
xlabel('Czas [s]'); ylabel('Poziom h_1 [cm]');
title('Poziom cieczy w pierwszym zbiorniku');
grid on;

subplot(3,1,3);
plot(t, h2, 'r', 'LineWidth', 1.5);
xlabel('Czas [s]'); ylabel('Poziom h_2 [cm]');
title('Poziom cieczy w drugim zbiorniku');
grid on;
print(fullfile(folderPath, 'Obiekt nieliniowy - zmiana zaklocenia.png'), '-dpng', '-r400');

% stabilizacja
T_end2 = 50000;
t2 = 0:dt:T_end2;
F_in_h2 = F1_0 * ones(size(t2));
FD_h2 = FD_0 * ones(size(t2));
F_in_h2(t2 > 1000) = 1.2 * F1_0;
F_in_h2_delay = [F1_0 * ones(1, N_delay), F_in_h2(1:end - N_delay)];

x = [h1_0; h2_0];
h1_h2 = zeros(size(t2));
h2_h2 = zeros(size(t2));

for k = 1:length(t2)
    h1_h2(k) = x(1);
    h2_h2(k) = x(2);
    dx = model_nieliniowy(t2(k), x, F_in_h2_delay(k), FD_h2(k));
    x = x + dx * dt;
end

figure('Name', 'Stabilizacja zmiennych układu', 'Color', 'w');
subplot(3,1,1);
plot(t2, F_in_h2, 'k--', 'LineWidth', 1.2); hold on;
plot(t2, F_in_h2_delay, 'r', 'LineWidth', 1.5);
xlabel('Czas [s]'); ylabel('F_{in} [cm^3/s]');
title('Sterowanie F_{in}');
legend('Wejście zadane','Wejście opóźnione');
grid on;

subplot(3,1,2);
plot(t2, h1_h2, 'b', 'LineWidth', 1.5);
xlabel('Czas [s]'); ylabel('Poziom h_1 [cm]');
title('Poziom cieczy w pierwszym zbiorniku');
grid on;

subplot(3,1,3);
plot(t2, h2_h2, 'r', 'LineWidth', 1.5);
xlabel('Czas [s]'); ylabel('Poziom h_2 [cm]');
title('Poziom cieczy w drugim zbiorniku');
grid on;
print(fullfile(folderPath, 'Obiekt nieliniowy - stabilizacja.png'), '-dpng', '-r400');

% linearyzacja 

disp('LINEARYZACJA MODELU W PUNKCIE PRACY');

T_end = 100000;
dt = 1;
t = 0:dt:T_end;

Fin = F1_0 * ones(size(t));
FD = FD_0 * ones(size(t));

Fin(t > 1000 & t <= T_end/2) = 1.1 * F1_0;
Fin(t > T_end/2) = 1.2 * F1_0;

N_delay = round(tau / dt);
Fin_delay = [F1_0 * ones(1, N_delay), Fin(1:end - N_delay)];

model_nieliniowy = @(t, x, Fin, FD) [
    (Fin + FD - alpha1*sqrt(x(1))) / (2*C1*x(1));
    (alpha1*sqrt(x(1)) - alpha2*sqrt(x(2))) / (3*C2*x(2)^2)
];

x = [h1_0; h2_0];
h1 = zeros(size(t));
h2 = zeros(size(t));

for k = 1:length(t)
    h1(k) = x(1);
    h2(k) = x(2);
    dx = model_nieliniowy(t(k), x, Fin_delay(k), FD(k));
    x = x + dx * dt;
end

% linearyzacja symboliczna
syms h1_sym h2_sym Fin_sym FD_sym real
syms alpha1_sym alpha2_sym C1_sym C2_sym real

f1 = (Fin_sym + FD_sym - alpha1_sym * sqrt(h1_sym)) / (2 * C1_sym * h1_sym);
f2 = (alpha1_sym * sqrt(h1_sym) - alpha2_sym * sqrt(h2_sym)) / (3 * C2_sym * h2_sym^2);

disp(' ');
disp('=== Nieliniowe równania różniczkowe ===');
fprintf('dh1/dt = %s\n', char(f1));
fprintf('dh2/dt = %s\n', char(f2));

A_sym = jacobian([f1; f2], [h1_sym, h2_sym]);
B_sym = jacobian([f1; f2], Fin_sym);

A = double(subs(A_sym, {h1_sym, h2_sym, Fin_sym, FD_sym, alpha1_sym, alpha2_sym, C1_sym, C2_sym}, ...
                {h1_0, h2_0, F1_0, FD_0, alpha1, alpha2, C1, C2}));
B = double(subs(B_sym, {h1_sym, h2_sym, Fin_sym, FD_sym, alpha1_sym, alpha2_sym, C1_sym, C2_sym}, ...
                {h1_0, h2_0, F1_0, FD_0, alpha1, alpha2, C1, C2}));
C_matrix = [0 1];
D = 0;


syms dh1_lin dh2_lin dFin real
f_lin = A_sym * [dh1_lin; dh2_lin] + B_sym * dFin;

disp(' ');
disp('=== Zlinearyzowany model w postaci równań (przybliżenie Taylora 1. rzędu) ===');
fprintf('d(dh1)/dt = %s\n', char(simplify(f_lin(1))));
fprintf('d(dh2)/dt = %s\n', char(simplify(f_lin(2))));


sys_lin = ss(A, B, C_matrix, D);

% Symulacja modelu liniowego dla 2 skoków
u_lin = Fin_delay' - F1_0;
[y_lin, ~, x_lin] = lsim(sys_lin, u_lin, t);
y_lin = y_lin + h2_0;

% wykresy
figure('Name', 'Porównanie modeli: nieliniowy vs liniowy (100000s)', 'Color', 'w');

subplot(2,1,1);
plot(t, Fin, 'k--', 'LineWidth', 1.5); hold on;
plot(t, Fin_delay, 'r', 'LineWidth', 1.5);
xlabel('Czas [s]');
ylabel('F_{in} [cm^3/s]');
title(['Sterowanie F_{in}(t) i opóźnione F_{in}(t - \tau), \tau = ' num2str(tau) ' s']);
legend('F_{in}(t)', 'F_{in}(t - \tau)', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t, h2, 'r', 'LineWidth', 1.5); hold on;
plot(t, y_lin, 'b--', 'LineWidth', 1.5);
xlabel('Czas [s]');
ylabel('Poziom h_2 [cm]');
title('Porównanie: model nieliniowy i zlinearyzowany (pełna dynamika)');
legend('Model nieliniowy', 'Model liniowy', 'Location', 'best');
grid on;

print(fullfile(folderPath, 'Porownanie_modeli_nieliniowy_liniowy_pelna_dynamika.png'), '-dpng', '-r400');

blad = h2 - y_lin';
fprintf('\nŚredni błąd aproksymacji: %.4f cm\n', mean(abs(blad)));
fprintf('Maksymalny błąd aproksymacji: %.4f cm\n', max(abs(blad)));


%Główne symulacje

T_end_comp = 50000;
dt_comp = 1;
t_comp = 0:dt_comp:T_end_comp;

skoki = [-10, 10, -20, 20];
nazwy_skokow = {'minus_10', 'plus_10', 'minus_20', 'plus_20'};

for i = 1:length(skoki)
    dF = skoki(i);
    
    Fin_comp = F1_0 * ones(size(t_comp));
    FD_comp = FD_0 * ones(size(t_comp));
    
    Fin_comp(t_comp > 1000) = F1_0 + dF;

    Fin_delay_comp = [F1_0 * ones(1, N_delay), Fin_comp(1:end - N_delay)];
    
    x_nl = [h1_0; h2_0];
    h1_nl = zeros(size(t_comp));
    h2_nl = zeros(size(t_comp));
    
    for k = 1:length(t_comp)
        h1_nl(k) = x_nl(1);
        h2_nl(k) = x_nl(2);
        dx = model_nieliniowy(t_comp(k), x_nl, Fin_delay_comp(k), FD_comp(k));
        x_nl = x_nl + dx * dt_comp;
    end
    
    u_lin_comp = Fin_delay_comp' - F1_0;
    [y_lin_comp, ~, x_lin_comp] = lsim(sys_lin, u_lin_comp, t_comp);
    h2_lin_comp = y_lin_comp' + h2_0;
    
    error_comp = h2_nl - h2_lin_comp;
    max_error = max(abs(error_comp));
    mean_error = mean(abs(error_comp));
    
    fprintf('Skok: %+d cm³/s | Błąd maks: %.4f cm | Błąd średni: %.4f cm\n', ...
            dF, max_error, mean_error);
    
    figure('Name', sprintf('Porównanie modeli - skok %+d', dF), 'Color', 'w', 'Position', [100 100 1000 700]);
    
    subplot(3,1,1);
    plot(t_comp, Fin_comp, 'k-', 'LineWidth', 1.5); hold on;
    plot(t_comp, Fin_delay_comp, 'r--', 'LineWidth', 1.2);
    ylabel('F_{in} [cm³/s]');
    title(sprintf('Sterowanie wejściowe - skok: %+d cm³/s', dF));
    legend('F_{in}(t)', 'F_{in}(t-τ) opóźnione', 'Location', 'best');
    grid on;
    xlim([0 50000]);
    
    subplot(3,1,2);
    plot(t_comp, h2_nl, 'b-', 'LineWidth', 2); hold on;
    plot(t_comp, h2_lin_comp, 'r--', 'LineWidth', 1.5);
    ylabel('h₂ [cm]');
    title('Odpowiedź wyjściowa - poziom h₂');
    legend('Model nieliniowy', 'Model liniowy', 'Location', 'best');
    grid on;
    xlim([0 50000]);
    
    subplot(3,1,3);
    plot(t_comp, error_comp, 'g-', 'LineWidth', 1.5);
    ylabel('Błąd [cm]');
    xlabel('Czas [s]');
    title('Błąd aproksymacji liniowej');
    grid on;
    xlim([0 50000]);
    print(fullfile(folderPath, sprintf('Porownanie_modeli_skok_%s.png', nazwy_skokow{i})), '-dpng', '-r400');
end