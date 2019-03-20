clear;
clc;
close all;
warning off;

N = 199;

sugeno_constant = readfis('lab_04_01.fis');
sugeno_linear   = readfis('lab_04_03.fis');

[Q, T, Qc, Qh] = calculate_plc(sugeno_constant, N);
[Ql, Tl, Qcl, Qhl] = calculate_plc(sugeno_linear, N);
 
figure();

%% Графики по константным выходным переменным

subplot(3,2,5);
plot(1:N+1,Qc,'b-');
hold on;
plot(1:N+1,Qh,'r-');
title('Qh, Qc');
xlabel('t,c');
ylabel('Q, %');
legend('Cold','Hot');
grid on;

subplot(3,2,1);
hold on;
plot(1:N,Q,'g-');
title('Q выходной');
xlabel('t,c');
ylabel('Q, %');
grid on;
hold on;
borders(N,50);

subplot(3,2,3)
plot(1:N,T,'b-');
title('Температура выходная');
xlabel('t,c');
ylabel('T, C');
grid on;
hold on;
borders(N,45);

%% Графики по линейным выходным переменным

subplot(3,2,6);
plot(1:N+1,Qcl,'b-');
hold on;
plot(1:N+1,Qhl,'r-');
title('Qh, Qc');
xlabel('t,c');
ylabel('Q, %');
legend('Cold','Hot');
grid on;

subplot(3,2,2);
hold on;
plot(1:N,Ql,'g-');
title('Q выходной');
xlabel('t,c');
ylabel('Q, %');
grid on;
hold on;
borders(N,50);

subplot(3,2,4)
plot(1:N,Tl,'b-');
title('Температура выходная');
xlabel('t,c');
ylabel('T, C');
grid on;
hold on;
borders(N,40);

%% Расчет переходных процессов
function [Q, T, Qc, Qh] = calculate_plc(fuzzy_control, N)

Qc = zeros(1,N);
Qh = zeros(1,N);

Q = zeros(1,N);
T = zeros(1,N);

Qc(1) = 20; %10^(-5);
Qh(1) = 20; %10^(-5);


dt = 1;

for i = 1:N
  
  t = i*dt;
  
  Tc = 10  + 10 * sin(pi*t/50 + 20);
  Th = 70 + 20 * cos(pi*t/60);

  T(i) = (Tc*Qc(i) + Th*Qh(i)) / (Qc(i) + Qh(i)); 
  Q(i) = (Qc(i) + Qh(i)) / 2;    

  input = [T(i), Q(i)];

  Tcurr = T(i);
  Qcurre = Q(i);
  
  dQhc = evalfis(input, fuzzy_control);

  Qc(i+1) = Qc(i) + dQhc(1);
  Qh(i+1) = Qh(i) + dQhc(2);
  
  Qccurr = Qc(i);
  Qhcurr = Qh(i);
    
end

end


%% Рисовалка трубки точности
function borders(x,y)
    plot([0 x],[0.90*y,0.90*y],'k--');
    plot([0 x],[1.10*y,1.10*y],'k--');
end

%%
%
%
%
%
%
%
%
%
%
%
%
%

