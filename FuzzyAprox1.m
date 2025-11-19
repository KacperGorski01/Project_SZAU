clear; clc;

Xend = 20;

f = @(x) sin(2*x)./exp(0.2*x);
x = linspace(0, Xend, 1000);

NumOfFuzzySets = 40;
Xp = linspace(0, Xend, NumOfFuzzySets);
DX = Xp(2) - Xp(1);

mf = cell(1, NumOfFuzzySets);
for i = 1:NumOfFuzzySets
    wsp = 0.1;
    % a
    if i == 1
        a = Xp(i);
    else
        a = Xp(i) - (0.5 + wsp)*DX;
    end
    % d
    if i == NumOfFuzzySets
        d = Xp(i);
    else
        d = Xp(i) + (0.5 + wsp)*DX;
    end
    % b
    if i == 1
        b = Xp(i);
    else
        b = Xp(i) - (0.5 - wsp)*DX;
    end
    % c
    if i == NumOfFuzzySets
        c = Xp(i);
    else
        c = Xp(i) + (0.5 - wsp)*DX;
    end

    mf{i} = @(z) trapmf(z, [a b c d]);
end

fLin = cell(1, NumOfFuzzySets);
for i = 1:NumOfFuzzySets
    eps = 1e-8;
    df = (f(Xp(i) + eps) - f(Xp(i) - eps)) / (2*eps);
    fLin{i} = @(x) df * (x - Xp(i)) + f(Xp(i));
end

fFuzzy = @(x) fuzzy_eval(x, mf, fLin, NumOfFuzzySets);
function y = fuzzy_eval(x, mf, fLin, NumOfFuzzySets)
    num = 0;
    den = 0;
    for i = 1:NumOfFuzzySets
        num = num + mf{i}(x) .* fLin{i}(x);
        den = den + mf{i}(x);
    end
    y = num ./ den;
end


% Rysowanie funkcji przynależności
figure(1);
hold on;
plot(x, f(x));
plot(x, fFuzzy(x));
for i = 1:NumOfFuzzySets
    plot(x, mf{i}(x));
end



