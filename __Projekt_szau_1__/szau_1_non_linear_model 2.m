clear all;
C1 = 0.75;      
C2 = 0.55;     
alpha1 = 20;    
alpha2 = 20;    

tau = 50;       
F1in = 52;      
FD = 11;       

% Warunki początkowe poziomów wody
h1_0 = 9.9225;      
h2_0 = 9.9225;  


czas_symulacji = [0 150]; 


[t, h] = ode45(@(t, h) zbiorniki(t, h, F1in, FD, C1, C2, alpha1, alpha2, tau), czas_symulacji, [h1_0, h2_0]);


figure('Position', [100, 50, 800, 400]);
subplot(2,1,1);
plot(t, h(:,1), 'b');
title('Symulacja poziomów w zbiornikach (F1 = 52, FD = 11, h1 = 9.9225, h2 = 9.9225)');
xlabel('Czas (s)');
ylabel('h1 (cm)');
legend('h1', 'h2');

subplot(2,1,2);
plot(t, h(:,2), 'r');
xlabel('Czas (s)');
ylabel('h2 (cm)');
legend('h2');

print(gcf, '1.4a.png', '-dpng', '-r300', '-opengl');

% Definicja funkcji do obliczania pochodnych
function dhdt = zbiorniki(t, h, F1in, FD, C1, C2, alpha1, alpha2, tau)
    % Zmienne stanu
    h1 = h(1);
    h2 = h(2);
    
    % Sygnał sterujący F1 z uwzględnieniem opóźnienia tau
    %F1 = F1in * heaviside(t - tau); % Przepływ F1 z uwzględnieniem opóźnienia tau
    F1 = F1in;
    
    % Przepływy
    F2 = alpha1 * sqrt(h1); 
    F3 = alpha2 * sqrt(h2);

    % Równania różniczkowe dla poziomów wody w zbiornikach
    dh1_dt = (F1 + FD - F2) / (2 * C1 * h1); 
    dh2_dt = (F2 - F3) / (3 * C2 * h2^2);  

    % Wektor pochodnych
    dhdt = [dh1_dt; dh2_dt];
end

