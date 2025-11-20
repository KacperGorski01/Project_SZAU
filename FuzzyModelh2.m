%% Porównanie modelu nieliniowego i zlinearizowanego
clear; clc;

Tend = 50e5;
h0 = [170; 100];
F1 = @(t) 200 + 300*(t >= 4e5) - 350*(t >= 40e5);
FD = @(t) 100 + 100*(t >= 30e5) - 100 *(t >= 45e5);

% Model rozmyty bardzo dobrze odwzorowuje wyjście h2. Stan h1 jest odwzorowywany zdecydowanie gorzej, ale ponieważ będziemy robić regulator DMC gdzie h2 jest 
% wyjściem układu, to właśnie na tym nam zależy, bo chcemy umieć robić dobrą predykcję wyjścia h2.


%% ---------------------------------------------------------------
% Model nieliniowy dynamiczny

fNonlinear = @(t,h) [ ...
    ( F1(t) + FD(t) - 23 * sqrt(h(1))) / (0.7 * h(1)) ; ...
    ( 23 * sqrt(h(1)) - 30 * sqrt(h(2))) / (1.35 * h(2)^2) ...
];


%% ---------------------------------------------------------------
%  Model zlinearyzowany

F1_point = 200;
Fd_point = 100;
h1_point = ( 1/23 * (F1_point + Fd_point) )^2;
h2_point = ( 1/30 * (F1_point + Fd_point) )^2;

A = zeros(2,2);
A(1,1) = -1.42857142857143*(F1_point + Fd_point - 23*sqrt(h1_point))/h1_point^2 - 16.4285714285714/h1_point^(3/2);
A(1,2) = 0;
A(2,1) = 8.51851851851852/(sqrt(h1_point)*h2_point^2);
A(2,2) = -1.48148148148148*(23*sqrt(h1_point) - 30*sqrt(h2_point))/h2_point^3 - 11.1111111111111/h2_point^(5/2);

B = zeros(2,2);
B(1,1) = 1.42857142857143/h1_point;
B(1,2) = 1.42857142857143/h1_point;
B(2,1) = 0;
B(2,2) = 0;

fLinear = @(t,h) A * (h - [h1_point; h2_point]) + B * [F1(t) - F1_point; FD(t) - Fd_point];


%% ---------------------------------------------------------------
% Model rozmyty z lokalnymi modelami zlinearizowanymi

% Zbiory rozmyte dla wyjścia h2
h2Range = [100, 500];

NumOfFuzzySets = 5;
h2p = linspace(h2Range(1), h2Range(2), NumOfFuzzySets);
Dh2p = h2p(2) - h2p(1);

mf = cell(1, NumOfFuzzySets);
sigma = 0.4*Dh2p;
for i = 1:NumOfFuzzySets
    c = h2p(i);
    mf{i} = @(z) gaussmf(z, [sigma c]);
end

figure(1);
x = linspace(h2Range(1), h2Range(2), 1000);
hold on;
title('Funkcje przynależności dla wyjścia h2');
xlabel('h2');
for i = 1:NumOfFuzzySets
    plot(x, mf{i}(x));
end

% Lokalne modele zlinearizowane
fLin = cell(1, NumOfFuzzySets);
for i = 1:NumOfFuzzySets
    h2p_ = h2p(i);
    Fdp_ = 100;
    F1p_ = sqrt(h2p_)*30 - Fdp_;
    h1p_ = ( 1/23 * (F1p_ + Fdp_) )^2;

    A = zeros(2,2);
    A(1,1) = -1.42857142857143*(F1p_ + Fdp_ - 23*sqrt(h1p_))/h1p_^2 - 16.4285714285714/h1p_^(3/2);
    A(1,2) = 0;
    A(2,1) = 8.51851851851852/(sqrt(h1p_)*h2p_^2);
    A(2,2) = -1.48148148148148*(23*sqrt(h1p_) - 30*sqrt(h2p_))/h2p_^3 - 11.1111111111111/h2p_^(5/2);

    B = zeros(2,2);
    B(1,1) = 1.42857142857143/h1p_;
    B(1,2) = 1.42857142857143/h1p_;
    B(2,1) = 0;
    B(2,2) = 0;

    fLin{i} = @(h,F) A * (h - [h1p_; h2p_]) + B * (F - [F1p_; Fdp_]);
end

% Model rozmyty
fFuzzy = @(t,h) fuzzy_eval(t, h, mf, fLin, NumOfFuzzySets, F1, FD);
function dh = fuzzy_eval(t, h, mf, fLin, NumOfFuzzySets, F1, FD)
    F1t = F1(t);
    FDt = FD(t);

    num = zeros(2,1);
    den = 0;
    for i = 1:NumOfFuzzySets
        mu = mf{i}(h(2));
        num = num + mu * fLin{i}(h, [F1t; FDt]);
        den = den + mu;
    end
    dh = num / den;    
end


%% ---------------------------------------------------------------
% Symulacja modeli nieliniowego i rozmytego

[tn, hn] = ode45(fNonlinear, [0, Tend], h0, odeset('RelTol',1e-6));
[tl, hl] = ode45(fLinear, [0, Tend], h0, odeset('RelTol',1e-6));
[tf, hf] = ode45(fFuzzy, [0, Tend], h0, odeset('RelTol',1e-6));

% Rysowanie wyników
figure(2);
subplot(2,1,1);
hold on;
plot(tn, hn(:,1), 'b', 'LineWidth', 1.5); 
plot(tf, hf(:,1), 'r--', 'LineWidth', 1.5);
plot(tl, hl(:,1), 'g-.', 'LineWidth', 1.5);
legend('Model nieliniowy', 'Model rozmyty', 'Model zlinearyzowany', 'Location', 'northeast');
grid on
grid minor
subplot(2,1,2);
hold on;
plot(tn, hn(:,2), 'b', 'LineWidth', 1.5); 
plot(tf, hf(:,2), 'r--', 'LineWidth', 1.5);
plot(tl, hl(:,2), 'g-.', 'LineWidth', 1.5);
legend('Model nieliniowy', 'Model rozmyty', 'Model zlinearyzowany', 'Location', 'northeast');
grid on;
grid minor










