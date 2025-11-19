clear; clc;

Xstart = 2;
Xend = 20;

f = @(x) 0.1*x.^2 - 0.2*(x-10).^2.*(x>10);
x = linspace(Xstart, Xend, 10000);

NumOfFuzzySets = 5;
Xp = linspace(Xstart, Xend, NumOfFuzzySets);
DX = Xp(2) - Xp(1);

mf = cell(1, NumOfFuzzySets);
sigma = 0.3*DX;
for i = 1:NumOfFuzzySets
    c = Xp(i);
    mf{i} = @(z) gaussmf(z, [sigma c]);
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



