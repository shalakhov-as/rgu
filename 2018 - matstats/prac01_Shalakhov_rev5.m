% 18.02.18 Исправлено:
% 1. T – вектор-столбец размерности [n x 1]
% 2. на выходе явно указаны возвращаемые параметры

function main()
clc;
close all;
tau = 60;
M = 120;
n = 10000;

% генерация одной выборки
T = generate_movement_time(10,tau,M,n);
ET = mean(T)
SKO = std(T)
DT = SKO^2

% построение гистограмм и плотностей распределения
figure;
hold on;

for k = 1:5
    % генерируем выборку для выбранного k
    T = generate_movement_time(k,tau,M,n);

    % строим гистограмму в левом столбике графиков
    subplot(5,2,(2*k - 1));
    hold on;
    t = title(['Гистограмма Т, ',num2str(k),' светофоров']);
    t.FontSize = 8;

    histogram(T);

    % строим распределения по хист_денсити
    % и по нормальному закону
    subplot(5,2,2*k);
    hold on;
    t = title(['Распределение вероятностей Т, ',num2str(k),' светофоров']);
    t.FontSize = 8;

    [p,x] = hist_density(T,20);
    plot(x,p);
    [p,x] = generate_normal_evaluation(T,k,tau,M,n);
    plot(x,p);
end

end

% -------------------- функции --------------------

% Генерирует n реализаций итогового Т, если на пути
% k светофоров с периодом tau и постоянная составляющая M
function T = generate_movement_time(k, tau, M, n)
    for i = 1:n
        T(i,1) = M + sum(rand(k,1) * tau);
    end
end

% Генерирует нормальное распределение итогового Т, если на пути
% k светофоров с периодом tau и постоянная составляющая M
function [p,x] = generate_normal_evaluation(T,k,tau,M,n)
    x = min(T):max(T);
    mu = M + tau * k / 2;
    sigma = sqrt( k * (tau^2) / 12 );
    p = normpdf(x,mu,sigma);
end
