clear; clc;

options = optimoptions('quadprog', "Algorithm","active-set"); options.Display = 'none'; 

% Zmienne modelu rozmytego

%liczba regulatorów
il_fun = 5;

h_min = 0;
h_max = 30;
h = (h_min:1:h_max)';


dumin = -inf;
dumax = inf;

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


%% Parametry regulatora

N = 400;
D = 1000;
lamb = 1;
Nu = 5;

Umax = 120;
Umin = 0;

lamb = lamb*ones(1,il_fun);

MP = zeros(N,D-1,il_fun);
M = zeros(N,Nu,il_fun);

for i = 1:il_fun
    s = generate_step_v2(hr0(i),false);
        
    for l=1:N
       for j=1:Nu
          if (l>=j)
             M(l,j,i)=s(l-j+1);
          end
       end
    end


    for l=1:N
       for j=1:D-1
          if l+j<=D
             MP(l,j,i)=s(l+j)-s(j);
          else
             MP(l,j,i)=s(D)-s(j);
          end    
       end
    end

end


%% Parametry modelu i symulacji

%Constants
C1 = 0.75;
C2 = 0.55;
ap1 = 20; %alfa_1
ap2 = 20; %alfa_2

%Punkt pracy obiektu
tau = 50;
h2_0 = 9.9225;
h1_0 = 9.9225;
F1 = 52;
FD = 11;

%Podstawowe punkty linearyzacji
F10 = 52;
FD0 = 11;

%Objetośc punktu pracy
v2_0 = h2_0^3 * C2;
v1_0 = h1_0^2 * C1;

t_sym = 8000; %czas symulacji
T = 1; %krok

ku = zeros(il_fun,D-1);
ke = zeros(1,il_fun);


%% Symulacja obiektu

%warunki_początkowe
kp = 500/T + 2;
kk = t_sym/T;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
h2_sl(1:kp) = h2_0;
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
Du = zeros(1,Nu)';
DUp = zeros(1,D-1)';
err_cur = 0;
err_sum = 0;
Yz = zeros(1,N)';
yk = zeros(1,N)';

A = [tril(ones(Nu));tril(ones(Nu))*-1]; %?
B = zeros(2*Nu,1); %?
ws = optimwarmstart(zeros(Nu,1),options);
%główne wykonanie programu
for k=kp:kk
    %symulacja obiektu
    v1(k) = real(v1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1*sqrt(h1(k-1))));
    v2(k) = real(v2(k-1) + T*(ap1*sqrt(abs(h1(k-1))) - ap2*(sqrt(abs(h2_sl(k-1))))));
    h1(k) = (v1(k) / C1)^(1/2);
    h2_sl(k) = (v2(k) / C2)^(1/3);
    
    %Liczenie błędu
    err_cur = yzad(k) - h2_sl(k);
    err_sum = err_sum + norm((yzad(k) - h2_sl(k)))^2;


    %Liczenie wartości przyrostu sterowania
    for i = 1:il_fun
        if i == 1
            w(i) = trapmf(h2_sl(k),[0 0 c(1)-nach/2 c(1)+ nach/2]);
        elseif i == il_fun
            w(i) = trapmf(h2_sl(k),[c(il_fun-1)-nach/2 c(il_fun-1)+nach/2 h_max h_max]);
        else
            w(i) = trapmf(h2_sl(k),[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]);
        end
    end
    
    %Sprawdzenie liczenia wag
    w_over_time(:,k) = w;

    Mr = zeros(N,Nu);
    MPr = zeros(N,D-1);
    for i = 1:il_fun
        Mr = Mr + w(i)*M(:,:,i)/sum(w);
        MPr = MPr + w(i)*MP(:,:,i)/sum(w);
    end
    lambr = w*lamb'/sum(w);

    Yz(1:end)=yzad(k);
    yk(1:end)=h2_sl(k);
    H = 2*(Mr'*Mr + lambr*eye(Nu,Nu));
    H = (H+H')/2;
    f = -2*Mr'*(Yz-yk-MPr*DUp);
    J = tril(ones(Nu,Nu));
    U = ones(Nu,1)*F1in(k-1);
    A_opt = [-J;J];
    B_opt = [Umin+U;Umax-U];
   test = quadprog(H,f,A_opt,B_opt,[],[],ones(Nu,1)*dumin, ones(Nu,1)*dumax, ws);
%     OPTIONS = optimoptions('fmincon','UseParallel',true, 'MaxFunctionEvaluations', 150);
%     Du = fmincon(@(Du)(Yz-yk-MPr*DUp-Mr*Du)'*(Yz-yk-MPr*DUp-Mr*Du)+lambr*Du'*Du,Du,A,B,[],[],ones(Nu,1)*-60,ones(Nu,1)*60);
   Du = test.X;
   holder = Du(1);
    for i = D-1:-1:2
        DUp(i) = DUp(i-1);
    end
    DUp(1) = holder;
    F1in(k) = F1in(k-1) + DUp(1);
    
%     figure
%     plot(h2, 'b')
%     hold on
%     plot(F1in,'g')
%     plot(yzad,"--r")
%     xlabel("k")
%     legend("Wyjście regulatora", "Sterowanie")
%     drawnow;

%     disp(k)
end

figure('Position', [100, 50, 800, 400]);
stairs(h2_sl, 'LineWidth', 2)
hold on;
stairs(yzad,"--", 'LineWidth', 2);
hold off;
xlabel('k'); ylabel("h");
legend("h_2","h_2_z_a_d")
title("FDMC")
print(gcf, 'h2.png', '-dpng', '-r300', '-opengl')

figure('Position', [100, 50, 800, 400]);
stairs(F1in)
legend("F_1_i_n")
xlabel('k'); ylabel("F_1_i_n");
legend("F1_in")
title("FDMC sterowanie")
print(gcf, 'F1.png', '-dpng', '-r300', '-opengl')

% display(err_sum)
save("SL_norm")
%%
function [s] = generate_step_v2(Y0, draw)
% Parametry programu
il_fun = 5;
C1 = 0.75;
C2 = 0.55;
ap1 = 20; %alfa_1
ap2 = 20; %alfa_2
 
%Punkt pracy
tau = 50;
h2_0 = 9.9225;
h1_0 = 9.9225;
v2_0 = h2_0^2 * C2;
F1 = 52;
FD = 11;

F10 = 52;
FD0 = 11;

t_sym = 8000; %czas symulacji
T = 1; %krok

% Wyliczanie sterowania co doprowadzi do Y0
s_h2_0 = Y0;
s_h1_0 = s_h2_0 * (ap2/ap1)^2;
s_FD0 = 11;
Wart_F1 = ap1*s_h1_0^0.5 - s_FD0;

%% Zmienne modelu rozmytego
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
m = (ap2/ap1)^2;
hr01 = hr0.*m;
vr1 = hr01.^2 * C1;
vr2 = hr0.^3 * C2;
Fr0 = ap1*hr01.^0.5-FD0;

%% Symulacja

if draw
    figure
    title("Odpowiedz skokowa")
    xlabel('k')
    ylabel('h_2')
    hold on
end

kp = tau/T + 2;
kk = t_sym/T;
k_s = 50;


h1 = h1_0 * ones(il_fun+1,kk);
h2_sl = h2_0 * ones(il_fun+1,kk);
v1 = h1_0^2 * C1 * ones(il_fun+1,kk);
v2 = h2_0^3 * C2 * ones(il_fun+1,kk);

    F1in(1:kk) = F1;
    FDc(1:kk) = FD;
    
    for k = kp:kk
        
        if k/T > k_s
            F1in(k) = Wart_F1;
        end

        for i = 1:il_fun
            %Rownania modelu
            
            v1(i,k) = v1(il_fun+1,k-1) + T*(F1in(k-1-(tau/T)) - Fr0(i) + FDc(k-1) - FD0 - (ap1/(2*(sqrt(hr01(i)))))*(h1(il_fun+1,k-1)-hr01(i)));
            v2(i,k) = v2(il_fun+1,k-1) + T*((ap1/(2*(sqrt(hr01(i)))))*(h1(il_fun+1,k-1)-hr01(i)) - (ap2/(2*(sqrt(hr0(i)))))*(h2_sl(il_fun+1,k-1)-hr0(i))); 
            h2_sl(i,k) = hr0(i) + (v2(i,k) - vr2(i))*(1/(3*(vr2(i)^2*C2)^(1/3)));
            h1(i,k) = hr01(i) + (v1(i,k) - vr1(i))*1/(2*sqrt(vr1(i)*C1));


            %Liczenie funkcji przynaleznosci
            if i == 1
                w(i) = trapmf(h2_sl(i,k),[0 0 c(1)-nach/2 c(1)+ nach/2]);
            elseif i == il_fun
                w(i) = trapmf(h2_sl(i,k),[c(il_fun-1)-nach/2 c(il_fun-1)+nach/2 h_max h_max]);
            else
                w(i) = trapmf(h2_sl(i,k),[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]);
            end
        end
        %Wyliczanie wyjscia modelu

        h2_sl(il_fun+1,k) = w * h2_sl(1:il_fun, k)/sum(w);
        h1(il_fun+1,k) = w * h1(1:il_fun, k)/sum(w);
        v2(il_fun+1,k) = w * v2(1:il_fun, k)/sum(w);
        v1(il_fun+1,k) = w * v1(1:il_fun, k)/sum(w);

    end

    s = h2_sl(6,k_s:end);
    s = s - s(1)*ones(1,length(s));
    s = s/(Wart_F1-F1);

% figure('Position', [100, 50, 800, 400]);
% stairs(h2, 'LineWidth', 2)
% hold on;
% stairs(yzad,"--", 'LineWidth', 2);
% hold off;
% xlabel('k'); ylabel("h");
% legend("h_2","h_2_z_a_d")
% title("FDMC")
% print(gcf, 'h2.png', '-dpng', '-r300', '-opengl')
% 
% figure('Position', [100, 50, 800, 400]);
% stairs(F1in)
% legend("F_1_i_n")
% xlabel('k'); ylabel("F_1_i_n");
% legend("F1_in")
% title("FDMC sterowanie")
% print(gcf, 'F1.png', '-dpng', '-r300', '-opengl')


end