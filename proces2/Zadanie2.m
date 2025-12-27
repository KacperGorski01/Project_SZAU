%% Projekt 2 - Sztuczna Inteligencja w Automatyce
% Zadanie 2
% 
% Punkt II - Modelowanie neuronowe przy użyciu sieci ELM
% Autorzy: Kacper Górski 314061, Antoni Marczuk 304436


%% 1. Wczytanie danych wygenerowanych w poprzednim kroku
clear all; close all;

if exist('dane_procesu.mat', 'file')
    load('dane_procesu.mat');
else
    error('Błąd: Nie znaleziono pliku dane_procesu.mat. Uruchom najpierw skrypt generujący dane.');
end

%% Przygotowanie danych w formacie macierzowym
% Tworzymy macierz wejść X, gdzie każdy wiersz to wektor:
% [u(k-4), u(k-5), y(k-1), y(k-2)]

d = 5; % Maksymalne opóźnienie (u(k-5))
N = length(u_train); % N = 2000
N_eff = N - d; % N_eff = 1995

X_train = zeros(N_eff, 4);
X_verify = zeros(N_eff, 4);
Y_train_target = y_train(d+1:end)';
Y_verify_target = y_verify(d+1:end)';

for k = (d + 1):N
    idx = k - d;
    % Wejścia sieci
    X_train(idx, :) = [u_train(k-4), u_train(k-5), y_train(k-1), y_train(k-2)];
    X_verify(idx, :) = [u_verify(k-4), u_verify(k-5), y_verify(k-1), y_verify(k-2)];
end

%% Główna pętla badania wpływu liczby neuronów K
K_values = 5:5:40; % K = 5, 10, 15, 20, 25, 30, 35, 40
liczba_powtorzen = 5; % Każda sieć uczona 5 razy

% Inicjalizacja tablic błędów
best_mse_train_step = zeros(size(K_values));
best_mse_verify_step = zeros(size(K_values));
best_mse_verify_rec = zeros(size(K_values));
inputSize = size(X_train, 2); 

% Struktura do przechowywania parametrów najlepszego wybranego modelu 
best_overall_model = struct('mse_rec', inf);

fprintf("--- Analiza modeli ELM ---\n");
fprintf("%-5s | %-10s | %-12s | %-12s | %-12s\n", 'K', 'Parameters', 'MSE_train', 'MSE_verify', 'MSE_rec_verify');
fprintf("----------------------------------------------------------------------\n");

% Dla każdego K:
for k_idx = 1:length(K_values)
    K = K_values(k_idx);
    current_best_rec_err = inf;
    best_p = 0; 
    rep_errors = zeros(1, liczba_powtorzen); 
    numParams = K*(inputSize+1) + (K+1); % Liczba parametrów
    params_list(k_idx) = numParams;
    
    % Pętla powtórzeń dla każdego K 
    for p = 1:liczba_powtorzen
        % 1. Losowanie wag warstwy ukrytej
        w10 = 2*(rand(K,1)-0.5);
        w1 = 2*(rand(K,inputSize)-0.5);
        
        % 2. Macierz wyjść warstwy ukrytej V
        V = tanh(w10 + w1 * X_train')'; 
        V_train = [ones(N_eff, 1), V];
        
        % 3. Obliczenie wag wyjściowych (lewe dzielenie)
        weightsW2 = V_train \ Y_train_target;
        w20 = weightsW2(1);
        w2 = weightsW2(2:end)';
        
        % 4. Modele w trybie bez rekurencji (jeden krok do przodu) 
        Ymod_train_step = model(w10, w1, w20, w2, X_train');
        Ymod_verify_step = model(w10, w1, w20, w2, X_verify');
        
        % 5. Model w trybie z rekurencją
        Ymod_verify_rec = zeros(1, N_eff);
        y_hist = [y_verify(d-1), y_verify(d)]; % y(k-2), y(k-1)
        
        for k = 1:N_eff
            x_rec = [u_verify(k+d-4); u_verify(k+d-5); y_hist(2); y_hist(1)];
            y_out = model(w10, w1, w20, w2, x_rec);
            Ymod_verify_rec(k) = y_out;
            y_hist = [y_hist(2), y_out];
        end
        
        mse_rec = mean((Y_verify_target' - Ymod_verify_rec).^2);
        rep_errors(p) = mse_rec; 
        
        % 6. Wybór najlepszego powtórzenia dla danego K 
        if mse_rec < current_best_rec_err
            current_best_rec_err = mse_rec;
            best_p = p;
            best_mse_verify_rec(k_idx) = mse_rec;
            best_mse_train_step(k_idx) = mean((Y_train_target' - Ymod_train_step).^2);
            best_mse_verify_step(k_idx) = mean((Y_verify_target' - Ymod_verify_step).^2);
            
            % Zapamiętanie globalnie najlepszego modelu
            if mse_rec < best_overall_model.mse_rec
                best_overall_model.mse_rec = mse_rec;
                best_overall_model.params = m_params(w10, w1, w20, w2, K);
            end
        end
    end
    
    % Wyświetlanie tabeli powtórzeń dla bieżącej liczby neuronów
    fprintf('\nWyniki dla K = %d (Parametry: %d)\n', K, numParams);
    fprintf('Nr próby | MSE Rekurencyjny (Weryfikacja)\n');
    fprintf('---------|-------------------------------\n');
    for p = 1:liczba_powtorzen
        if p == best_p
            fprintf('   %d     | %e [WYBRANY]\n', p, rep_errors(p));
        else
            fprintf('   %d     | %e\n', p, rep_errors(p));
        end
    end
end

%% Wyświetlanie tabeli dla każdej liczby neuronów
fprintf("\n\n=== PODSUMOWANIE ZBIORCZE ===\n");
fprintf("%-5s | %-8s | %-12s | %-12s | %-12s\n", 'K', 'Params', 'MSE_Train', 'MSE_Verify', 'MSE_Rec_Ver');
fprintf("--------------------------------------------------------------------\n");
for k_idx = 1:length(K_values)
    fprintf("%-5d | %-8d | %-12.2e | %-12.2e | %-12.2e\n", ...
        K_values(k_idx), params_list(k_idx), best_mse_train_step(k_idx), ...
        best_mse_verify_step(k_idx), best_mse_verify_rec(k_idx));
end
fprintf("--------------------------------------------------------------------\n");

[~, idx_opt] = min(best_mse_verify_rec);
m = best_overall_model.params; 

fprintf("----------------------------------------------------------------------\n");
fprintf("Wybrany model ostateczny: K = %d\n", m.K);
fprintf("Błąd rekurencyjny: %e\n\n", best_mse_verify_rec(idx_opt));

%% Wykresy dla wybranego modelu
datasets = {struct('u', u_train, 'y', y_train, 'name', 'Uczący'), ...
            struct('u', u_verify, 'y', y_verify, 'name', 'Weryfikujący')};

for i = 1:2
    [y_rec, y_step] = final_sim(m, datasets{i}.u, datasets{i}.y, d, N);
    figure;
    subplot(2,1,1);
    plot(datasets{i}.y(d+1:end), 'k'); hold on; plot(y_step, 'r--');
    title(['Zbiór ', datasets{i}.name, ' (K=', num2str(m.K), ') - Bez rekurencji']); grid on;
    legend('Dane', 'Model');
    
    subplot(2,1,2);
    plot(datasets{i}.y(d+1:end), 'k'); hold on; plot(y_rec, 'b--');
    title(['Zbiór ', datasets{i}.name, ' - Tryb rekurencyjny']); grid on;
    xlabel('k'); legend('Dane', 'Model');
end

%% funkcje pomocnicze
function Ymod = model(w10, w1, w20, w2, X)
    Ymod = w20+w2*tanh(w10+w1*X);
end

function p = m_params(w10, w1, w20, w2, K)
    p = struct('w10', w10, 'w1', w1, 'w20', w20, 'w2', w2, 'K', K);
end

function [y_rec, y_step] = final_sim(m, u, y, d, N)
    N_eff = N - d;
    y_rec = zeros(1, N_eff);
    y_step = zeros(1, N_eff);
    y_hist = [y(d-1), y(d)];
    for k = 1:N_eff
        % Krok w przód
        x_step = [u(k+d-4); u(k+d-5); y(k+d-1); y(k+d-2)];
        y_step(k) = model(m.w10, m.w1, m.w20, m.w2, x_step);
        % Rekurencja
        x_rec = [u(k+d-4); u(k+d-5); y_hist(2); y_hist(1)];
        y_out = model(m.w10, m.w1, m.w20, m.w2, x_rec);
        y_rec(k) = y_out;
        y_hist = [y_hist(2), y_out];
    end
end