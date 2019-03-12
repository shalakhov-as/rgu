clear;
clc;
close all;
warning on;

t = 0;
tmax = 600;

% 40C is comfortable temperature
% 50% is comfortable flow

Qc = zeros(1,tmax); % cold flow
Qh = zeros(1,tmax); % hot flow

Q = zeros(1,tmax); % result flow
T = zeros(1,tmax); % result temperature

Qc(1) = 10^(-5); % initial cold flow
Qh(1) = 10^(-5); % initial hot flow

Q(1) = 0; % initial result flow
T(1) = 0; % initial result temperature

dt = 1;
i = 2;

fuzzy_control = readfis('lab_02.fis');

while(t <= tmax)
  
  Tc = 10  + 10 * sin(pi*t/50 + 20);
  Th = 70 + 20 * cos(pi*t/60);

  T(i) = (Tc*Qc(i-1) + Th*Qh(i-1)) / (Qc(i-1) + Qh(i-1)); 
  Q(i) = (Qc(i-1) + Qh(i-1)) / 2;    

  input = [T(i), Q(i)];

  Tcurr = T(i);
  Qcurre = Q(i);
  
  dQhc = evalfis(input, fuzzy_control);

  Qc(i) = Qc(i-1) + dQhc(2);
  if (Qc(i) > 100)
    error('Значение вышло за пределы диапазона');
    %Qc(i) = 100;
  end
  if (Qc(i) < 0)
    error('Значение вышло за пределы диапазона');
    %Qc(i) = 0;
  end
  
  Qh(i) = Qh(i-1) + dQhc(1);
  if (Qh(i) > 100)
    error('Значение вышло за пределы диапазона');
    %Qh(i) = 100;
  end
  if (Qh(i) < 0)
    error('Значение вышло за пределы диапазона');
    %Qh(i) = 0;
  end
  
  Qccurr = Qc(i);
  Qhcurr = Qh(i);

  i = i+1;
  t = t+dt;
    
end

figure;
hold on;
title('Result flow over time');
plot(1:t+1, Q, 'LineWidth',2);
borders(t, 50);
grid on;

figure;
hold on;
title('Result temperature over time');
plot(1:t+1, T, 'LineWidth',2);
borders(t,40);
grid on;

figure;
hold on;
title('Crane flow over time');
plot(1:t+1, Qc, 'b', 'LineWidth',2);
plot(1:t+1, Qh, 'r', 'LineWidth',2);
grid on;


function borders(x,y)
    plot([0 x],[0.9*y,0.9*y],'k--');
    plot([0 x],[1.1*y,1.1*y],'k--');
end
