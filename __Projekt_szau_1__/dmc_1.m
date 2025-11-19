% --- Parametry DMC ---
Nu = 5;            % Horyzont sterowania
N = 400;            % Horyzont predykcji
Du = 1000;          
lambda = 1;         % Lambda

ks = 1500;          % Początek sterowania
kp = 1400;          % Punkt przejściowy
kk = 8000;          % Liczba kroków symulacji
T = 1;              % Krok dyskretyzacji

% --- Punkt linearyzacji ---
C1 = 0.75;          % Stała dla zbiornika 1
C2 = 0.55;          % Stała dla zbiornika 2
alpha1 = 20;        % Współczynnik przepływu z 1 do 2
alpha2 = 20;        % Współczynnik odpływu
F1 = 52;

% --- Symulacja odpowiedzi skokowej ---
F1_0 = 52;          % Skok wejściowy
FD_0 = 11;          % Zakłócenie
h1_0 = 9.9225;      % Początkowy poziom h1
h2_0 = 9.9225;      % Początkowy poziom h2
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

% --- Regulatora DMC ---
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

% --- Symulacja systemu z regulacją DMC ---
h1_dmc = h1_0;
h2_dmc = h2_0;
v1_dmc = v1_0;
v2_dmc = v2_0;
h2_dmc_values = zeros(1, numel(step_time));
delta_u = zeros(1, numel(step_time));
h2_zad(1:1000)=9.9225; 
h2_zad(1000:3000)=14;
h2_zad(3000:8000)=26;

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
stairs(step_time, h2_zad, 'b');
hold on;
stairs(step_time, h2_dmc_values(1:k), 'r');
xlabel('Czas [s]');
ylabel('h2 [cm]');
title('Regulator DMC analityczny');
legend('Wartość zadana h2', 'Wyjście h2');
grid on;
print(gcf, '1.16.png', '-dpng', '-r300', '-opengl');

%% 



set(0,'DefaultStairLineWidth',1);
Umax = 120;
Umin = 0;

Nu = 5;
N = 400;
D = 1000;
lambda = 1;

C1 = 0.75;
C2 = 0.55;
alpha1 = 20;
alpha2 = 20;
F1 = 52; 
FD = 11; 
tau = 50; 
h1_eq = 9.9225;
h2_eq = 9.9225;

%liczba regulatorów
il_fun = 5;
lambda = lambda*ones(1,il_fun);


%Punkt pracy obiektu
h2_0 = h2_eq;
h1_0 = h1_eq;


%Podstawowe punkty linearyzacji
F10 = 52;
FD0 = 11;

%Objetośc punktu pracy
v2_0 = h2_0^3 * C2;
v1_0 = h1_0^2 * C1;

t_sym = 8000; %czas symulacji
T = 1; %krok


h_min = 0;
h_max = 30;
h = (h_min:1:h_max)';

nach = 1; %nachylenie funkcji 

d = (h_max-h_min)/il_fun; %szerokości funkcji przynależnośći
c = h_min+d:d:h_max-d; %punkty przegięcia

%Wybranie punktu linearyzacji
hr0 = ones(1,il_fun);
hr0(1) = d/2;
hr0(il_fun) = min((h_max+c(il_fun-1))/2+1, h_max);
if il_fun > 2
    hr0(2:il_fun-1) = (c(2:il_fun-1)+c(1:il_fun-2))./2;
end


ku = zeros(il_fun,D-1);
ke = zeros(1,il_fun);


for r = 1:il_fun
    s = generate_step_fuzzy_DMC_1(hr0(r), C1, C2, alpha1, alpha2, tau);
    k_s(:,:,r) = s;
    M=zeros(N,Nu);
        for i=1:N
           for j=1:Nu
              if (i>=j)
                 M(i,j)=s(i-j+1);
              end
           end
        end
        

    MP=zeros(N,D-1);
        for i=1:N
           for j=1:D-1
              if i+j<=D
                 MP(i,j)=s(i+j)-s(j);
              else
                 MP(i,j)=s(D)-s(j);
              end      
           end
        end
    
    K = ((M'*M + lambda(r) * eye(Nu))^(-1))* M';
    ku(r,:) = K(1,:)*MP;
    ke(r) = sum(K(1,:));
end

%warunki_początkowe
kp = tau/T + 2;
ks = max(19,500+100); %chwila skoku wartosci zadania
kk = t_sym/T;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
h2(1:kp) = h2_0;
h1(1:kp) = h1_0;
F1in(1:T:kp) = F1;
FD = 11;
FDc(1:T:t_sym/T) = FD;

%Skok wartosci zadanej:
yzad(1:1000)=9.9225; 
yzad(1000:3000)=14;
yzad(3000:8000)=26;


error = 0;
w = zeros(1,il_fun);
Du = zeros(il_fun,1);
DUp = zeros(1,D-1);
Y = zeros(N,1);
err_cur = 0;

DUfin = 0;

%główne wykonanie programu
for k=kp:kk
    %symulacja obiektu
    v1(k) = real(v1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - alpha1*sqrt(h1(k-1))));
    v2(k) = real(v2(k-1) + T*(alpha1*sqrt(abs(h1(k-1))) - alpha2*(sqrt(abs(h2(k-1))))));
    h1(k) = (v1(k) / C1)^(1/2);
    h2(k) = (v2(k) / C2)^(1/3);

  
    %Liczenie błędu
    err_cur = yzad(k) - h2(k);

    for i = D-1:-1:2
      DUp(i) = DUp(i-1);
    end

    DUp(1) = DUfin;

    %Liczenie wartości przyrostu sterowania
    for i = 1:il_fun

        Du(i) = ke(i)*err_cur-ku(i,:)*DUp';

        if i == 1
            w(i) = trapmf(h2(k),[0 0 c(1)-nach/2 c(1)+ nach/2]);
        elseif i == il_fun
            w(i) = trapmf(h2(k),[c(il_fun-1)-nach/2 c(il_fun-1)+nach/2 h_max h_max]);
        else
            w(i) = trapmf(h2(k),[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]);
        end
    end
    
    %Sprawdzenie liczenia wag
%     w_over_time(:,k) = w;

    %Ogranieczenia przyrostu sterowania
    DUfin = w * Du / sum(w);
    

    F1in(k) = F1in(k-1) + DUfin;

    %Ograniczenia sterowania
    if F1in(k) > Umax
        F1in(k) = Umax;
    elseif F1in(k) < Umin
        F1in(k) = Umin;
    end
    
    DUfin = F1in(k) - F1in(k-1); 

end


iteracja = 0:1:kk-1;
%Plot wyjście
figure('Position', [100, 50, 800, 400]);
stairs(iteracja, h2, 'LineWidth', 2)
hold on;
stairs(iteracja, yzad,"--", 'LineWidth', 2);
hold off;
xlabel('k'); ylabel("h");
legend("h_2","h_2_z_a_d")
title("FDMC")
print(gcf, 'h2.png', '-dpng', '-r300', '-opengl')
odpowiedz_dmc_fuzyy = h2;

%Plot sterowanie
figure('Position', [100, 50, 800, 400]);
stairs(iteracja, F1in)
legend("F_1_i_n")
xlabel('k'); ylabel("F_1_i_n");
legend("F1_in")
title("FDMC sterowanie")
print(gcf, 'F1.png', '-dpng', '-r300', '-opengl')

%% 
figure('Position', [100, 50, 800, 400]);
stairs(iteracja, h2, 'LineWidth', 2)
hold on;
stairs(iteracja, yzad,"g--", 'LineWidth', 1);
hold on;
stairs(step_time(2:k), h2_dmc_values(2:k), 'r', 'LineWidth', 2);
hold on;
stairs(h2_sl,"--", 'LineWidth', 2)
xlabel('k'); ylabel("h");
legend("h_2 FDMC","h_2 zad", "h_2 DMC","h2_sl")
title("Porównanie FDMC i DMC")
print(gcf, 'porownaniemodeli.png', '-dpng', '-r300', '-opengl')
