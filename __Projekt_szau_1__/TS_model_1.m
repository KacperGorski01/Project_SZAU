function TS_model_1(nach, ilosc_modeli, h1_eq, h2_eq, C1, C2, alpha1, alpha2, tau, F1, FD)

set(0,'DefaultStairLineWidth',1);


%Punkt pracy

h2_0 = h2_eq;
h1_0 = h1_eq;
v2_0 = h2_0^3 * C2;
v1_0 = h1_0^2 * C1;

F10 = 52;
FD0 = 11;

t_sym = 10000; %czas symulacji
T = 1 ; %krok


h_min = 0;
h_max = 30;
nach = 1;
ilosc_modeli = 2;

d = (h_max-h_min)/ilosc_modeli; %szerokości funkcji przynależnośći
c = h_min+d:d:h_max-d; %punkty przegięcia

%Wybranie punktu linearyzacji
hr0 = ones(1,ilosc_modeli);
hr0(1) = d/2;
hr0(ilosc_modeli) = min((h_max+c(ilosc_modeli-1))/2+1, h_max);
    if ilosc_modeli > 2
        hr0(2:ilosc_modeli-1) = (c(2:ilosc_modeli-1)+c(1:ilosc_modeli-2))./2;
    end
m = (alpha2/alpha1)^2;
hr01 = hr0.*m;
vr2 = hr0.^3 * C2;
vr1 = hr01.^2 * C1;
Fr0 = alpha1*hr01.^0.5-FD0;


    figure
    hold on
    %Plotter funkcji przynaleznosci
    for i = 1:ilosc_modeli
        if i == 1
            plot(trapmf(h_min:1:h_max,[0 0 c(1)-nach/2 c(1)+ nach/2]), 'LineWidth', 2);
        elseif i == ilosc_modeli
            plot(trapmf(h_min:1:h_max,[c(ilosc_modeli-1)-nach/2 c(ilosc_modeli-1)+nach/2 h_max h_max]), 'LineWidth', 2);
        else
            plot(trapmf(h_min:1:h_max,[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]), 'LineWidth', 2);
        end
    end
    xlim([0 30])
    xlabel("h_2"); ylabel("Funkcja przynależności");
    title(sprintf("Funkcja przynaleznosci dla %i zbiorów rozmytych i nachylenia %f",ilosc_modeli, nach))
    
legends = {};


    figure
    hold on
    % title("Porownanie obiektów")
    xlabel('k')
    ylabel('h_2')


kp = tau/T + 2;
kk = t_sym/T;

colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']; 

for P = 32:10:72


    h1 = h1_0 * ones(ilosc_modeli+3,kk);
    h2 = h2_0 * ones(ilosc_modeli+3,kk);
    v1 = h1_0^2 * C1 * ones(ilosc_modeli+3,kk);
    v2 = h2_0^3 * C2 * ones(ilosc_modeli+3,kk);

    F1in(1:kk) = F1;
    FDc(1:kk) = FD;
    
    for k = kp:kk
        
        if k/T > 50
            F1in(k) = P;
        end

        % Model nieliniowy
        v1(ilosc_modeli+2,k) = v1(ilosc_modeli+2,k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - alpha1*(sqrt(h1(ilosc_modeli+2,k-1))));
        v2(ilosc_modeli+2,k) = v2(ilosc_modeli+2,k-1) + T*(alpha1*sqrt(h1(ilosc_modeli+2,k-1)) - alpha2*(sqrt(h2(ilosc_modeli+2,k-1)))); 
        h1(ilosc_modeli+2,k) = sqrt(v1(ilosc_modeli+2,k)/C1);
        h2(ilosc_modeli+2,k) = (v2(ilosc_modeli+2,k)/C2)^(1/3);

        % Model liniowy
        v1(ilosc_modeli+3,k) = v1(ilosc_modeli+3,k-1) + T*(F1in(k-1-(tau/T)) - F10 + FDc(k-1) - FD0 - (alpha1/(2*(sqrt(h1_0))))*(h1(ilosc_modeli+3,k-1)-h1_0));
        v2(ilosc_modeli+3,k) = v2(ilosc_modeli+3,k-1) + T*((alpha1/(2*(sqrt(h1_0))))*(h1(ilosc_modeli+3,k-1)-h1_0) - (alpha2/(2*(sqrt(h2_0))))*(h2(ilosc_modeli+3,k-1)-h2_0)); 
        h2(ilosc_modeli+3,k) = h2_0 + 1/(3*(v2_0^2*C2)^(1/3))*(v2(ilosc_modeli+3,k) - v2_0);
        h1(ilosc_modeli+3,k) = h1_0 + 1/(2*sqrt(v1(i,k)*C1))*(v1(ilosc_modeli+3,k) - v1_0);

        for i = 1:ilosc_modeli
            %Rownania modelu
       
            v1(i,k) = v1(ilosc_modeli+1,k-1) + T*(F1in(k-1-(tau/T)) - Fr0(i) + FDc(k-1) - FD0 - (alpha1/(2*(sqrt(hr01(i)))))*(h1(ilosc_modeli+1,k-1)-hr01(i)));
            v2(i,k) = v2(ilosc_modeli+1,k-1) + T*((alpha1/(2*(sqrt(hr01(i)))))*(h1(ilosc_modeli+1,k-1)-hr01(i)) - (alpha2/(2*(sqrt(hr0(i)))))*(h2(ilosc_modeli+1,k-1)-hr0(i))); 
            h2(i,k) = hr0(i) + (v2(i,k) - vr2(i))*(1/(3*(vr2(i)^2*C2)^(1/3)));
            h1(i,k) = hr01(i) + (v1(i,k) - vr1(i))*1/(2*sqrt(v1(i,k)*C1));

            %Liczenie funkcji przynaleznosci
            if i == 1

                w(i) = trapmf(h2(i,k),[0 0 c(1)-nach/2 c(1)+ nach/2]);
            elseif i == ilosc_modeli
                w(i) = trapmf(h2(i,k),[c(ilosc_modeli-1)-nach/2 c(ilosc_modeli-1)+nach/2 h_max h_max]);
            else
                w(i) = trapmf(h2(i,k),[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]);
            end
        end
        %Wyliczanie wyjscia modelu

        h2(ilosc_modeli+1,k) = w * h2(1:ilosc_modeli, k)/sum(w);
        h1(ilosc_modeli+1,k) = w * h1(1:ilosc_modeli, k)/sum(w);
        v2(ilosc_modeli+1,k) = w * v2(1:ilosc_modeli, k)/sum(w);
        v1(ilosc_modeli+1,k) = w * v1(1:ilosc_modeli, k)/sum(w);

    end

% Przedstawienie wyników modeli
    colorIndex = mod(find(P == [49 70.2 78 85.8 117]) - 1, length(colors)) + 1;
    stairs((1:k),h2(ilosc_modeli+2,:),['--', colors(colorIndex)])
    stairs((1:k),h2(ilosc_modeli+3,:),['-.', colors(colorIndex)])
    stairs((1:k),h2(ilosc_modeli+1,:),['-', colors(colorIndex)])
    
    % Dodanie etykiety dla bieżącej wartości P
    legends{end+1} = ['F1in = ', num2str(P), ' (Model nieliniowy)'];
    legends{end+1} = ['F1in = ', num2str(P), ' (Model zlinearyzowany)'];
    legends{end+1} = ['F1in = ', num2str(P), ' (Model rozmyty)'];

    clear v1 v2 h1 h2
end

    legend(legends, 'Location', 'best', 'Orientation', 'vertical');

    title("Porownanie modeli rozmytych, nachylenie = " + nach + ", liczba modeli = " + ilosc_modeli)

clear stairs legends w F1in d c 

end
