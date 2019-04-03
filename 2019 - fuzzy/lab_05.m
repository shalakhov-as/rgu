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
steps = 500;

disp('Расчет нечёткого регулирования с синусоидальными возмущениями');
fuzzy_controller();
disp('Расчет нечёткого регулирования с случайными возмущениями');
fuzzy_agress();
disp('Расчет ПИД-регулирования с случайными возмущениями');
pid_controller();


function fuzzy_controller()

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
  g2(1)=0.5;

  P0 = 0.35;
  P(1) = P0;
  E(1) = P0-SP;
  t = dt;

  fuzzy_control = readfis('lab_05.fis');

  for i=2:steps

    dg2 = evalfis([P(i-1), g2(i-1)], fuzzy_control);
    g2(i) = g2(i-1) + dg2;

    t = t+dt;

    P1(i) = 0.5*(0.05*sin(t)+1);
    P2(i) = 0.3*(0.05*sin(t+10*dt)+1);

    Ps(i) = findPs(P1(i),P2(i),g2(i-1));
    ro = P(i-1)/(z*R*T);
    a(i) = get_a(P1(i),P2(i),Ps(i),g2(i-1),tau,P0,S,ro,ksi);
    P(i) = (P(i-1)-Ps(i))*exp(-a(i)*dt)+Ps(i);
    E(i) = (P(i)-SP);

  end

  figure(1);
  subplot(3,2,1);
  plot(0:dt:t-dt,P,'b-'), xlabel('time'), ylabel('P');
  title('fuzzyreg');
  axis([0 steps*dt 0.3 0.6]);
  grid on;
  hold on
  borders(t-dt,SP);

  figure(1);
  subplot(3,2,2);
  plot(0:dt:t-dt,g2,'r-'), xlabel('time'), ylabel('gamma');
  axis([0 steps*dt 0 1]);
  title('fuzzyreg');
  grid on;
    

end

function fuzzy_agress()
  
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
  g2(1)=0.5;

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
  P1(1)=0.5;
  P2(2)=0.3;

  for i=2:steps
    
      dg2 = evalfis([P(i-1),g2(i-1)],fuzzy_control);
      g2(i) = g2(i-1)+dg2;
      t = t+dt;        
      P1(i) = 0.5 + 0.07*P1(i-1) + e1(i);
      P2(i) = 0.3 - 0.05*P2(i-1) + e2(i);
      Ps(i) = findPs(P1(i),P2(i),g2(i-1));
      ro = P(i-1)/(z*R*T);
      a(i) = get_a(P1(i),P2(i),Ps(i),g2(i-1),tau,P0,S,ro,ksi);
      P(i) = (P(i-1)-Ps(i))*exp(-a(i)*dt)+Ps(i);
      E(i) = (P(i)-SP);
      
  end

  figure(1);
  subplot(3,2,3);
  plot(0:dt:t-dt,P,'r-'), xlabel('hh'), ylabel('hh');
  title('fuzzyreg, autoregress')
  axis([0 steps*dt 0.3 0.6]);
  grid on;
  hold on;
  
  borders(t-dt,SP);

  figure(1);
  subplot(3,2,4);
  plot(0:dt:t-dt,g2), xlabel('h'), ylabel('h');
  axis([0 steps*dt 0 1]);
  title('fuzzyreg, autoregress');
  grid on;

end

function pid_controller()
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
  
  t = 0;
  g2(1) = 0.5;

  P0 = 0.4;
  P(1) = P0;
  E(1) = P0-SP;
  t = dt;

  I = 0;
  kc = 5;
  ki = 1;
  
  for i=2:steps
    
    if i<=2
      [dg2,I]= PID(E(i-1),0,kc,ki,I,dt);
    else
      [dg2,I]= PID(E(i-1),E(i-2),kc,ki,I,dt);
    end
    
    g2(i) = dg2;
    t = t+dt;
    P1(i) = 0.5*(0.05*sin(t)+1);
    P2(i) = 0.3*(0.05*sin(t+10*dt)+1);
    Ps(i) = findPs(P1(i),P2(i),g2(i-1));
    ro = P(i-1)/(z*R*T);
    a(i) = get_a(P1(i),P2(i),Ps(i),g2(i-1),tau,P0,S,ro,ksi);
    P(i) = (P(i-1)-Ps(i))*exp(-a(i)*dt)+Ps(i);
    E(i) = (P(i)-SP);
    
  end

  figure(1);
  subplot(3,2,5);
  plot(0:dt:t-dt,P,'r-'), xlabel('hh'), ylabel('hh'), title('hh');
  axis([0 steps*dt 0.3 0.6]);
  grid on;
  hold on;
  borders(t-dt,SP);

  figure(1);
  subplot(3,2,6);
  plot(0:dt:t-dt,g2,'b-'), xlabel('hh'), ylabel('hh'), title('hh');
  axis([0 steps*dt 0 1]);
  grid on;
  
end

function Ps = findPs(P1,P2,g2)

  p=[-1 (P1-g2^2*P2) +g2^2*P2^2];
  c=roots(p);
  if c(1)==0||c(2)==0
      Ps=c(1);   
  else
      pos=find(c>0);
      Ps=c(pos);
  end
    
end

function a = get_a(P1,P2,Ps,g2,tau,P0,S,ro,ksi)

  M1 = 1*S*power(2*ro*P0*(P1-P0)/ksi,0.5);
  M2 = g2*S*power(2*ro*P2*(P0-P2)/ksi,0.5);
  dM = M1-M2;
  a = dM/(tau*(Ps-P0));
    
end

function [y,I] = PID(ec,em1,kc,ki,I,dt)

    D = (ec-em1);
    kd = 0;
    
    if ki~=0
      
        Imax = (1-(kc*ec+kd*D))/ki;
        Imin = -(kc*ec+kd*D)/ki;
        I = I + ec*dt;
        I = min(max(I,Imin), Imax);
        
    end

    OP = kc*ec+ki*I+kd*D/dt;
    
    if OP > 1
        OP = 1;
    elseif OP < 0
        OP = 0;
    end
    y = OP;
end

function borders(x,y)
  figure(1);
  plot([0 x],[0.95*y,0.95*y],'r');
  plot([0 x],[1.05*y,1.05*y],'r');
end