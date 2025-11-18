clear all;
C1 = 0.75;      
C2 = 0.55;      
alpha1 = 20;    
alpha2 = 20;    

tau = 50;      
F1in = 52;      
FD = 11;       

% Warunki początkowe poziomów cieczy
h1_0 = 7;      
h2_0 = 9.9225;  

% Warunki początkowe dla objętości V1 i V2
V1_0 = C1 * h1_0^2;      
V2_0 = C2 * h2_0^3;      

czas_symulacji = [0 150]; 


[t_h, h] = ode45(@(t, h) zbiorniki_h(t, h, F1in, FD, C1, C2, alpha1, alpha2, tau), czas_symulacji, [h1_0, h2_0]);
[t_V, V] = ode45(@(t, V) zbiorniki_V(t, V, F1in, FD, C1, C2, alpha1, alpha2, tau), czas_symulacji, [V1_0, V2_0]);

% Przeliczenie objętości na poziomy w zbiornikach
h1_V = sqrt(V(:, 1) / C1);
h2_V = nthroot(V(:, 2) / C2, 3);

% Rysowanie wyników symulacji
figure('Position', [100, 50, 800, 400]);

subplot(2, 1, 1);
plot(t_h, h(:, 1), 'b', t_V, h1_V, 'g--');
title('Poziom w zbiorniku 1 (h1)');
xlabel('Czas (s)');
ylabel('h1 (cm)');
legend('h1 - zmienne stanu wysokości', 'h1 - zmienne stanu objętości');

subplot(2, 1, 2);
plot(t_h, h(:, 2), 'r', t_V, h2_V, 'g--');
title('Poziom w zbiorniku 2 (h2)');
xlabel('Czas (s)');
ylabel('h2 (cm)');
legend('h2 - zmienne stanu wysokości', 'h2 - zmienne stanu objętości');

print(gcf, '1.4.png', '-dpng', '-r300', '-opengl');

% Definicja funkcji do obliczania pochodnych dla poziomów
function dhdt = zbiorniki_h(t, h, F1in, FD, C1, C2, alpha1, alpha2, tau)
    % Zmienne stanu
    h1 = h(1);
    h2 = h(2);
    
    
    F1 = F1in;
    
    % Przepływy
    F2 = alpha1 * sqrt(h1); 
    F3 = alpha2 * sqrt(h2); 

    % Równania różniczkowe dla poziomów
    dh1_dt = (F1 + FD - F2) / (2 * C1 * h1);
    dh2_dt = (F2 - F3) / (3 * C2 * h2^2);

    % Wektor pochodnych
    dhdt = [dh1_dt; dh2_dt];
end

% Definicja funkcji do obliczania pochodnych dla objętości
function dVdt = zbiorniki_V(t, V, F1, FD, C1, C2, alpha1, alpha2, tau)
    % Zmienne stanu (objętości)
    V1 = V(1);
    V2 = V(2);
    
    % Przeliczenie objętości na poziomy
    h1 = sqrt(V1 / C1);
    h2 = nthroot(V2 / C2, 3);
    
    % Przepływy
    F2 = alpha1 * sqrt(h1); 
    F3 = alpha2 * sqrt(h2); 
    
    % Równania różniczkowe dla objętości
    dV1_dt = F1 + FD - F2;
    dV2_dt = F2 - F3;
    
    % Wektor pochodnych
    dVdt = [dV1_dt; dV2_dt];
end
