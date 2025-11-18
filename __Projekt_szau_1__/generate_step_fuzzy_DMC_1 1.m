function [s] = generate_step_fuzzy_DMC_1(Y0, C1, C2, alpha1, alpha2, tau)


il_fun = 5;
set(0,'DefaultStairLineWidth',1);

h2_0 = 9.925;
h1_0 = 9.925;
v2_0 = h2_0^2 * C2;
F1 = 52;
FD = 11;

F10 = 52;
FD0 = 11;

t_sym = 10000; %czas symulacji
T = 1; %krok

s_h2_0 = Y0;
s_h1_0 = s_h2_0 * (alpha2/alpha1)^2;
s_v1_0 = s_h1_0^2 * C1;
s_v2_0 = s_h2_0^3 * C2;
s_FD0 = 11;
Wart_F1 = alpha1*s_h1_0^0.5 - s_FD0;

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
m = (alpha2/alpha1)^2;
hr01 = hr0.*m;
vr1 = hr01.^2 * C1;
vr2 = hr0.^3 * C2;
Fr0 = alpha1*hr01.^0.5-FD0;

kp = tau/T + 2;
kk = t_sym/T;
k_s = 150;


h1 = h1_0 * ones(il_fun+1,kk);
h2 = h2_0 * ones(il_fun+1,kk);
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
            
            v1(i,k) = v1(il_fun+1,k-1) + T*(F1in(k-1-(tau/T)) - Fr0(i) + FDc(k-1) - FD0 - (alpha1/(2*(sqrt(hr01(i)))))*(h1(il_fun+1,k-1)-hr01(i)));
            v2(i,k) = v2(il_fun+1,k-1) + T*((alpha1/(2*(sqrt(hr01(i)))))*(h1(il_fun+1,k-1)-hr01(i)) - (alpha2/(2*(sqrt(hr0(i)))))*(h2(il_fun+1,k-1)-hr0(i))); 
            h2(i,k) = hr0(i) + (v2(i,k) - vr2(i))*(1/(3*(vr2(i)^2*C2)^(1/3)));
            h1(i,k) = hr01(i) + (v1(i,k) - vr1(i))*1/(2*sqrt(vr1(i)*C1));


            %Liczenie funkcji przynaleznosci
            if i == 1
                w(i) = trapmf(h2(i,k),[0 0 c(1)-nach/2 c(1)+ nach/2]);
            elseif i == il_fun
                w(i) = trapmf(h2(i,k),[c(il_fun-1)-nach/2 c(il_fun-1)+nach/2 h_max h_max]);
            else
                w(i) = trapmf(h2(i,k),[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]);
            end
        end
        %Wyliczanie wyjscia modelu

        h2(il_fun+1,k) = w * h2(1:il_fun, k)/sum(w);
        h1(il_fun+1,k) = w * h1(1:il_fun, k)/sum(w);
        v2(il_fun+1,k) = w * v2(1:il_fun, k)/sum(w);
        v1(il_fun+1,k) = w * v1(1:il_fun, k)/sum(w);

    end

    s = h2(6,k_s:end);
    s = s - s(1)*ones(1,length(s));
    s = s/(Wart_F1-F1);

end