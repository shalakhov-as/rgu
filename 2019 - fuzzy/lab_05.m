clc;
clear;
close all;

global S;
global ksi;
global z;
global R;
global T;
global V;
global SP;
global dt;
global tau;

global kp;
global ki;

global steps;

S = 0.3;
ksi = 15000;
z = 0.96;
R = 500;
T = 300;
V = 100;
SP = 0.45;
dt = 10;
tau = V/(z*R*T);
steps = 400;

kp = 3;
ki = 1;

disp('Расчет нечёткого регулирования с син. возмущениями');
fuzzy_reg_sin();
disp('Расчет нечёткого регулирования с автокорр. возмущениями');
fuzzy_reg_autocorr();
disp('Расчет ПИД-регулирования с син. возмущениями');
pid_reg_sin();
disp('Расчет ПИД-регулирования с авторег. возмущениями');
pid_reg_autocorr();


function fuzzy_reg_sin()

  global S;
  global ksi;
  global z;
  global R;
  global T;
  global V;
  global SP;
  global dt;
  global tau;
  global steps;

  t=0;
  gamma2(1)=0.5;

  P0 = 0.35;
  P(1) = P0;
  E(1) = P0-SP;
  t = dt;

  fuzzy_control = readfis('lab_05.fis');

  for i=2:steps

    d_gamma2 = evalfis([P(i-1), gamma2(i-1)], fuzzy_control);
    gamma2(i) = gamma2(i-1) + d_gamma2;

    t = t+dt;

    P1(i) = 0.5*(0.05*sin(t)+1);
    P2(i) = 0.3*(0.05*sin(t+10*dt)+1);

    Ps(i) = calculate_pressure_stat(P1(i),P2(i),gamma2(i-1));
    ro = P(i-1)/(z*R*T);
    a(i) = get_a(P1(i),P2(i),Ps(i),gamma2(i-1),tau,P0,S,ro,ksi);
    P(i) = (P(i-1)-Ps(i))*exp(-a(i)*dt)+Ps(i);
    E(i) = (P(i)-SP);

  end

  figure(1);
  subplot(4,2,1);
  plot(0:dt:t-dt,P, 'LineWidth', 2), xlabel('t,c'), ylabel('P, Мпа');
  title('fuzzy, син. возмущения');
  axis([0 steps*dt 0.3 0.6]);
  grid on;
  hold on;
  borders(t-dt,SP);

  figure(1);
  subplot(4,2,2);
  plot(0:dt:t-dt,gamma2, 'LineWidth', 2), xlabel('t,c'), ylabel('gamma 2, %');
  axis([0 steps*dt 0 1]);
  title('fuzzy, син. возмущения');
  grid on;
    
end

function fuzzy_reg_autocorr()
  
  global S;
  global ksi;
  global z;
  global R;
  global T;
  global V;
  global SP;
  global dt;
  global tau;
  global steps;
  
  t=0;
  gamma2(1)=0.5;

  P0=0.4;
  P(1) = P0;
  E(1) = P0-SP;
  t = dt;

  fuzzy_control = readfis('lab_05.fis');

  sigma1 = 0.01;
  sigma2 = 0.012;
  mu = 0.01; 
  e1 = normrnd(mu,sigma1,[1 steps]);
  e2 = normrnd(mu,sigma2,[1 steps]);
  P1(1) = 0.5;
  P2(1) = 0.3;

  for i=2:steps
    
      d_gamma2 = evalfis([P(i-1),gamma2(i-1)],fuzzy_control);
      gamma2(i) = gamma2(i-1)+d_gamma2;
      t = t+dt;        
      P1(i) = 0.5 + 0.07*P1(i-1) + e1(i);
      P2(i) = 0.3 - 0.05*P2(i-1) + e2(i);
      Ps(i) = calculate_pressure_stat(P1(i),P2(i),gamma2(i-1));
      ro = P(i-1)/(z*R*T);
      a(i) = get_a(P1(i),P2(i),Ps(i),gamma2(i-1),tau,P0,S,ro,ksi);
      P(i) = (P(i-1)-Ps(i))*exp(-a(i)*dt)+Ps(i);
      E(i) = (P(i)-SP);
      
  end

  figure(1);
  subplot(4,2,3);
  plot(0:dt:t-dt,P, 'LineWidth', 2), xlabel('t,c'), ylabel('P, Мпа');
  title('fuzzy, авторег. возмущения');
  axis([0 steps*dt 0.3 0.6]);
  grid on;
  hold on;
  
  borders(t-dt,SP);

  figure(1);
  subplot(4,2,4);
  plot(0:dt:t-dt,gamma2, 'LineWidth', 2), xlabel('t,c'), ylabel('gamma 2, %');
  axis([0 steps*dt 0 1]);
  title('fuzzy, авторег. возмущения');
  grid on;

end

function pid_reg_sin()
  global S;
  global ksi;
  global z;
  global R;
  global T;
  global V;
  global SP;
  global dt;
  global tau;
  global steps;
  global kp;
  global ki;
  
  t = 0;
  gamma2(1) = 0.5;

  P0 = 0.4;
  P(1) = P0;
  E(1) = P0-SP;
  t = dt;

  I = 0;
  
  for i=2:steps
    
    if i<=2
      [OP,I]= calculate_pid_control(E(i-1),0,kp,ki,I,dt);
    else
      [OP,I]= calculate_pid_control(E(i-1),E(i-2),kp,ki,I,dt);
    end
    
    gamma2(i) = OP;
    t = t+dt;
    P1(i) = 0.5*(0.05*sin(t)+1);
    P2(i) = 0.3*(0.05*sin(t+10*dt)+1);
    Ps(i) = calculate_pressure_stat(P1(i),P2(i),gamma2(i-1));
    ro = P(i-1)/(z*R*T);
    a(i) = get_a(P1(i),P2(i),Ps(i),gamma2(i-1),tau,P0,S,ro,ksi);
    P(i) = (P(i-1)-Ps(i))*exp(-a(i)*dt)+Ps(i);
    E(i) = (P(i)-SP);
    
  end

  figure(1);
  subplot(4,2,5);
  plot(0:dt:t-dt,P, 'LineWidth', 2), xlabel('t,c'), ylabel('P, Мпа');
  title('ПИД, син. возмущения');
  axis([0 steps*dt 0.3 0.6]);
  grid on;
  hold on;
  borders(t-dt,SP);

  figure(1);
  subplot(4,2,6);
  plot(0:dt:t-dt,gamma2, 'LineWidth', 2), xlabel('t,c'), ylabel('gamma 2, %');
  title('ПИД, син. возмущения');
  axis([0 steps*dt 0 1]);
  grid on;
  
end

function pid_reg_autocorr()
  global S;
  global ksi;
  global z;
  global R;
  global T;
  global V;
  global SP;
  global dt;
  global tau;
  global steps;
  global kp;
  global ki;
  
  t = 0;
  gamma2(1) = 0.5;

  P0 = 0.4;
  P(1) = P0;
  E(1) = P0-SP;
  t = dt;

  I = 0;

  sigma1 = 0.01;
  sigma2 = 0.012;
  mu = 0.01; 
  e1 = normrnd(mu,sigma1,[1 steps]);
  e2 = normrnd(mu,sigma2,[1 steps]);
  P1(1) = 0.5;
  P2(1) = 0.3;
  
  for i=2:steps
    
    if i<=2
      [OP,I]= calculate_pid_control(E(i-1),0,kp,ki,I,dt);
    else
      [OP,I]= calculate_pid_control(E(i-1),E(i-2),kp,ki,I,dt);
    end
    
    gamma2(i) = OP;
    t = t+dt;
    P1(i) = 0.5 + 0.07*P1(i-1) + e1(i);
    P2(i) = 0.3 - 0.05*P2(i-1) + e2(i);
    Ps(i) = calculate_pressure_stat(P1(i),P2(i),gamma2(i-1));
    ro = P(i-1)/(z*R*T);
    a(i) = get_a(P1(i),P2(i),Ps(i),gamma2(i-1),tau,P0,S,ro,ksi);
    P(i) = (P(i-1)-Ps(i))*exp(-a(i)*dt)+Ps(i);
    E(i) = (P(i)-SP);
    
  end

  figure(1);
  subplot(4,2,7);
  plot(0:dt:t-dt,P, 'LineWidth', 2), xlabel('t,c'), ylabel('P, Мпа');
  title('ПИД, авторег. возмущения');
  axis([0 steps*dt 0.3 0.6]);
  grid on;
  hold on;
  borders(t-dt,SP);

  figure(1);
  subplot(4,2,8);
  plot(0:dt:t-dt,gamma2, 'LineWidth', 2), xlabel('t,c'), ylabel('gamma 2, %');
  title('ПИД, авторег. возмущения');
  axis([0 steps*dt 0 1]);
  grid on;
  
end

function Ps = calculate_pressure_stat(P1,P2,gamma2)

  p = [-1 (P1-gamma2^2*P2) +gamma2^2*P2^2];
  c = roots(p);
  if c(1)==0 || c(2)==0
      Ps = c(1);   
  else
      i = find(c>0);
      Ps = c(i);
  end
    
end

function a = get_a(P1,P2,Ps,gamma2,tau,P0,S,ro,ksi)

  gamma1 = 1;

  M1 = gamma1 * S * (2*ro*P0*(P1-P0)/ksi) ^ 0.5;
  M2 = gamma2 * S * (2*ro*P2*(P0-P2)/ksi) ^ 0.5;
  dM = M1-M2;
  a = dM / (tau*(Ps-P0));
    
end

function [OP,I] = calculate_pid_control(e_current,e_prev,kp,ki,I,dt)

    D = (e_current-e_prev);
    kd = 0;
    
    if ki ~= 0
      
        Imax = (1-(kp*e_current+kd*D))/ki;
        Imin = -(kp*e_current+kd*D)/ki;
        I = I + e_current*dt;
        I = min(max(I,Imin), Imax);
        
    end

    OP = kp*e_current + ki*I + kd*D/dt;
    
    if OP > 1
        OP = 1;
    elseif OP < 0
        OP = 0;
    end

end

function borders(x,y)
  figure(1);
  plot([0 x],[0.95*y,0.95*y],'r--');
  plot([0 x],[1.05*y,1.05*y],'r--');
end