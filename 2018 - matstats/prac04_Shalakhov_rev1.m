% TODO рассчитать характеристики
% средневыборочных по формулам из лекции

function main()
  clc;
  close all;
  clear all;

  mu = 10;
  sko = 5;
  D = sko*sko;

  n = 10000;
  N = 10000;

  % 24 vectors of 1,000 elements each
  t = generate_randvals(n, N, mu, sko);
  % single vector of 1,000 elements
  t1 = t(:,1);

  % vector of means, 24 elements
  tMean = mean(t)';
  % expected value of means is the same
  % as of values themselves
  muMean = mu;
  % dispersion of means equals D/n
  skoMean = sqrt(D/n);

  % visualize probability density functions
  % of values and means of samples
  build_characteristics(t1, mu, sko, 'Random values');
  build_characteristics(tMean, muMean, skoMean, 'Means of random values');

  % build pdf of values in means' pdfs window
  [x,p] = hist_density(tMean, 20);
  pn = normpdf(x,mu,sko);
  plot(x,pn,'--','LineWidth', 2,'DisplayName','normpdf, given values');
end

function build_characteristics(t, mu, sko, figName)
  [x,p] = hist_density(t, 20);
  pn = normpdf(x,mu,sko);

  figure('NumberTitle', 'off', 'Name', figName);
  subplot(2,1,1);
  hold on;
  title('Histogram');
  histogram(t,50,'DisplayName','histogram of T');

  subplot(2,1,2);
  hold on;
  title('Probability density finctions')
  plot(x,pn,'--','LineWidth', 2,'DisplayName','normpdf');
  plot(x,p,'DisplayName','HistDensity');

  legend('show');
end

function x = generate_randvals(n, N, E, sko)
  x = randn(n,N)*sko + E;
end

function [centers, density] = hist_density(x, bin_count)
  dx = (max(x) - min(x)) / bin_count;
  [counts, centers] = hist(x, bin_count);
  density = (counts / length(x)) / dx;
end
