function main()
  clc;
  close all;
  clear all;

  n = 4;
  N = 10000;
  D = 4;
  SKO = sqrt(D);

  T = randn(n,N).*SKO;

  xn = mean(T);
  % z-statistics
  z = xn./sqrt(D/n);
  % t-statistics
  t = xn./(std(T,1)/sqrt(n));

  % remove values that are too great or
  % too small to build more detailed
  % hist_density in target range of values
  % i.e. remove all "100" or "-200" values
  ind = find(abs(t) < 3);
  t = t(ind);

  % plotting

  % task 4-5
  [xnx, xnp] = hist_density(xn,30);
  xnpt = normpdf(xnx, 0, sqrt(D/n));

  figure;

  subplot(4,1,1);
  hold on;
  legend;
  title('Mean(x) distribution')
  plot(xnx,xnp,'DisplayName','Evaluation');
  plot(xnx,xnpt,'DisplayName','Theoretical');
  xlim([-3 3]);

   % task 6
  [zx, zp] = hist_density(z,30);
  zpt = normpdf(zx, 0, 1);

  subplot(4,1,2);
  hold on;
  legend;
  title('Z-stats distribution');
  plot(zx,zp,'DisplayName','Evaluation');
  plot(zx,zpt,'DisplayName','Theoretical');
  xlim([-3 3]);

  % task 7
  [tx, tp] = hist_density(t,30);
  tpt = tpdf(tx, n-1);

  subplot(4,1,3);
  hold on;
  legend;
  title('T-stats distribution')
  plot(tx,tp,'DisplayName','Evaluation');
  plot(tx,tpt,'DisplayName','Theoretical');
  xlim([-3 3]);

  % task 8
  % Find n value which makes T-statistics look like
  % Z-statistics

  % Answer:
  % with n = 20 T-stats look almost like Z-stats

  n = 20;
  T = randn(n,N).*SKO;
  xn = mean(T);
  z = xn./sqrt(D/n);
  t = xn./(std(T,1)/sqrt(n));

  [tx, tp] = hist_density(t,30);
  [zx, zp] = hist_density(z,30);

  subplot(4,1,4);
  hold on;
  legend;
  title('T-stats vs Z-stats comparizon, n=20')
  plot(tx,tp,'DisplayName','T');
  plot(zx,zp,'DisplayName','Z');
  xlim([-3 3]);

end

function [centers, density] = hist_density(x, bin_count)
  dx = (max(x) - min(x)) / bin_count;
  [counts, centers] = hist(x, bin_count);
  density = (counts / length(x)) / dx;
end
