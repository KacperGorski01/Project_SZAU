clear all;


C1 = 0.75;      
C2 = 0.55;      
alpha1 = 20;    
alpha2 = 20;    

tau = 50;      
F1in = 52;      
FD = 40;       

% Warunki początkowe
h1_0 = 0.9225;      
h2_0 = 9.9225;  

h1_1 = h1_0;
h2_2 = h2_0;

% Czas symulacji
czas_symulacji = [0 300]; 


[t_h, h] = ode45(@(t, h) zbiorniki_h(t, h, F1in, FD, C1, C2, alpha1, alpha2, tau), czas_symulacji, [h1_0, h2_0]);
[t_h2, h2] = ode45(@(t, h) zbiorniki_lin(t, h, F1in, FD, C1, C2, alpha1, alpha2, tau), czas_symulacji, [h1_1, h2_2]);


figure('Position', [100, 50, 800, 400]);

subplot(2, 1, 1);
plot(t_h, h(:, 1), 'b', t_h2, h2(:,1), 'g--');
title('Symulacja poziomów w zbiornikach (F1 = 52, FD = 40, h1 = 0.9225, h2 = 9.9225)');
xlabel('Czas (s)');
ylabel('h1 (cm)');
legend('Model nieliniowy', 'Model liniowy', 'Location', 'best');  

subplot(2, 1, 2);
plot(t_h, h(:, 2), 'r', t_h2, h2(:,2), 'g--');
xlabel('Czas (s)');
ylabel('h2 (cm)');
legend('Model nieliniowy', 'Model liniowy', 'Location', 'best');  


%print(gcf, '1.15.png', '-dpng', '-r300', '-opengl');

function dhdt = zbiorniki_h(t, h, F1in, FD, C1, C2, alpha1, alpha2, tau)
    % Zmienne stanu
    h1 = h(1);
    h2 = h(2);
    
    F1 = F1in;
    
    % Przepływy
    F2 = alpha1 * sqrt(h1); 
    F3 = alpha2 * sqrt(h2); 

    % Równania dla zmiany poziomów
    dh1_dt = (F1 + FD - F2) / (2 * C1 * h1);
    dh2_dt = (F2 - F3) / (3 * C2 * h2^2);

    % Wektor pochodnych
    dhdt = [dh1_dt; dh2_dt];
end


function dhdt2 = zbiorniki_lin(t, h, F1in, FD, C1, C2, alpha1, alpha2, tau)
    % Zmienne stanu
    h1 = h(1);
    h2 = h(2);
    
    % Punkt pracy
    h1_ss = 9.9225;
    h2_ss = 9.9225;

    F1 = F1in;
    
    
    dh1_dt = (F1 + FD - alpha1 * (sqrt(h1_ss) + (1 / (2 * sqrt(h1_ss))) * (h1 - h1_ss))) / (2 * C1 * h1_ss);
    dh2_dt = ((alpha1 * (sqrt(h1_ss) + (1 / (2 * sqrt(h1_ss))) * (h1 - h1_ss))) - ...
              (alpha2 * (sqrt(h2_ss) + (1 / (2 * sqrt(h2_ss))) * (h2 - h2_ss)))) / ...
             (3 * C2 * (h2_ss^2 + 2 * h2_ss * (h2 - h2_ss)));

    % Wektor pochodnych
    dhdt2 = [dh1_dt; dh2_dt];
end

