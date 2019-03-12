function main()
  clc;
  clear all;
  close all;
  warning off;
  n = 60;
  [t, p] = aperiodic(n);
  [p1, p2] = generate_sensor_data(60, 0.1, 0.1, 0.2, 0.1);

  figure;

  subplot(2,1,1);
  hold on;
  legend on;
  title('P(t), kg/cm2');
  plot(t,p,'--','LineWidth', 2,'DisplayName','P');
  plot(t,p1,'DisplayName','P1');
  plot(t,p2,'DisplayName','P2');

  subplot(2,1,2);
  title('P1 - P2');
  hold on;
  legend on;
  plot(t,p1 - p2,'DisplayName','diff');

  % E(dP) = systematic(P1) - systematic(P2)
  % D(dP) = D(P1) + D(P2)
end

function [P1, P2] = generate_sensor_data(n, systematic1, dispersion1, systematic2, dispersion2)
  [t,p] = aperiodic(n);
  noise = randn(2,n);
  P1 = p + noise(1,:) * sqrt(dispersion1) + systematic1;
  P2 = p + noise(2,:) * sqrt(dispersion2) + systematic2;
end

function [t,p] = aperiodic(n)
  t = 1:n;
  tau = (1/n)*3;
  p = (exp(t.*(-1).*tau).*(-1) + 1).*20 + 10;
end
