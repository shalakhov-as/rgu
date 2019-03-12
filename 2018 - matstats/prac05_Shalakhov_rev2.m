% rev2
% randchi2: "transpose" removed
% dispersion_evaluation: simplified with std

function main()

  clc;
  clear all;
  close all;

  % part 1

  k = 5;
  ch = randchi2(k,1000000);
  [x,p] = hist_density(ch, 1000);
  pt = chi2pdf(x,k);

  figure;
  hold on;
  legend;
  title('Chi Square distrubutions comparison');
  plot(x,p,'DisplayName','Custom');
  plot(x,pt,'DisplayName','Built-in');
  xlim([0 20]);

  % part 2
  n = 5;
  DT = 3;
  SKO = sqrt(DT);

  % populating samples
  T = randn(n,100000).*SKO;

  % evaluate dispersion
  [Du, Ds] = dispersion_evaluation(T);
  [du,pu] = hist_density(Du, 1000);
  [ds,ps] = hist_density(Ds, 1000);

  figure;
  hold on;
  legend;
  title('Dispersion evaluations distribution comparison');
  plot(du,pu,'DisplayName','Unshifted');
  plot(ds,ps,'DisplayName','Shifted');

  xlim([0 10]);

  % part 2.4

  Dc = (Du./DT).*(n-1);
  Ch = randchi2(n-1, 100000);
  [dc,pc] = hist_density(Dc, 1000);
  [ch,ph] = hist_density(Ch, 1000);
  figure;
  hold on;
  legend;
  title('Evaluation of (n-1)S2/DT')
  plot(dc,pc,'DisplayName','(n-1) S2 / DT');
  plot(ch,ph,'DisplayName','chi2(n-1)');
  xlim([0 30]);

end

% generate sample of n elements
% of chi-square distribution
function x = randchi2(k,n)
  x = sum(randn(n,k).^2, 2);
end

% evaluates dispersions of samples in x
function [Du, Ds] = dispersion_evaluation(x)
  Ds = std(x,1).^2;
  Du = std(x,0).^2;
end

function [centers, density] = hist_density(x, bin_count)
  dx = (max(x) - min(x)) / bin_count;
  [counts, centers] = hist(x, bin_count);
  density = (counts / length(x)) / dx;
end
