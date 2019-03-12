function main()

  clc;
  close all;
  clear all;

  part3();
  part4(8);
  part4(16);
  part5();

end

function part3()

  % IMPORTANT NOTE!
  % P1 and P2 is defined as matrices, not as vectors!
  % Functions were made with this in mind
  %
  % Modification is approved by V.V.Yuzhanin

  D1 = 3;
  D2 = 3;

  n = 6;
  N = 1000000;

  a = 0.05;

  % E-statistics analisys
  figure;
  subplot(3,1,1);
  title('E-stats distributions');
  legend;
  xlim([-4 8]);
  ylim([0 0.5]);
  hold on;

  x = -5:0.1:5;
  p = normpdf(x, 0, sqrt((D1+D2)/n));
  plot(x, p, 'DisplayName','Theoretical E, H0');
  p = normpdf(x, 2, sqrt((D1+D2)/n));
  plot(x, p, 'DisplayName','Theoretical E, H1');

  sys11 = 2;
  sys21 = 0;
  [P11,P21] = generate_sensor_data(n, N, sys11, D1, sys21, D2);
  est = e_statistic(P11, P21);
  [x, p] = hist_density(est,50);
  plot(x, p, 'DisplayName','Practical E, H1');

  [P10,P20] = generate_sensor_data(n, N, 0, D1, 0, D2);
  est = e_statistic(P10, P20);
  [x, p] = hist_density(est,50);
  plot(x, p, 'DisplayName','Practical E, H0');

  [h, border_low, border_high] = test_equal_means_e(est, D1, D2, n, a);
  plot(border_low, 0,'or', border_high, 0, 'or');

  % Z-statistics analisys
  subplot(3,1,2);
  legend;
  title('Z-stats distributions');
  xlim([-4 8]);
  ylim([0 0.5]);
  hold on;

  x = -5:0.1:5;
  p = normpdf(x, 0, 1);
  plot(x, p, 'DisplayName','Theoretical Z, H0');
  x = x + (sys11-sys21)/sqrt((D1+D2)/n);
  plot(x, p, 'DisplayName','Theoretical Z, H1');

  zst = z_statistic(P11, P21, D1, D2);
  [x, p] = hist_density(zst,200);
  plot(x, p, 'DisplayName','Practical Z, H1');

  zst = z_statistic(P10, P20, D1, D2);
  [x, p] = hist_density(zst,200);
  plot(x, p, 'DisplayName','Practical Z, H0');

  [h, border_low, border_high] = test_equal_means_z(zst, a);
  plot(border_low, 0,'or', border_high, 0, 'or');

  % T-statistics analisys
  subplot(3,1,3);
  legend;
  title('T-stats distributions');
  xlim([-4 8]);
  ylim([0 0.5]);
  hold on;

  tst1 = t_statistic(P11, P21);
  tst0 = t_statistic(P10, P20);

  x = -5:0.1:5;
  p = tpdf(x, n-1);
  plot(x, p, 'DisplayName','Theoretical T, H0');
  %x = x + (sys11 - sys21)/mean(std(P1-P2,1)/sqrt(size(P1,1)));
  %plot(x, p, 'DisplayName','Theoretical T, H1');

  [x, p] = hist_density(tst1,50);
  plot(x, p, 'DisplayName','Practical T, H1');
  [x, p] = hist_density(tst0,50);
  plot(x, p, 'DisplayName','Practical T, H0');

  [h, border_low, border_high] = test_equal_means_t(tst0, n, a);
  plot(border_low, 0,'or', border_high, 0, 'or');

end

function part4(n)

  N = 10000;
  D1 = 3;
  D2 = 3;

  disp('PART 4, n = ');
  disp(n);

  % 4.1 - Calculate error of first kind
  %       in case of main hypothesis
  [P1,P2] = generate_sensor_data(n, N, 0, D1, 0, D2);

  est = e_statistic(P1, P2);
  [h, border_low, border_high] = test_equal_means_e(est, D1, D2, n, 0.05);
  disp('(Part 4.1) Error of 1st kind (Estat):');
  disp(1-h);

  zst = z_statistic(P1, P2, D1, D2);
  [h, border_low, border_high] = test_equal_means_z(zst, 0.05);
  disp('(Part 4.3) Error of 1st kind (Zstat):');
  disp(1-h);

  % 4.1 - Calculate error of second kind
  %       in case of secondary hypothesis
  [P1,P2] = generate_sensor_data(n, N, 2, D1, 0, D2);

  est = e_statistic(P1, P2);
  [h, border_low, border_high] = test_equal_means_e(est, D1, D2, n, 0.05);
  disp('(Part 4.2) Error of 2nd kind (Estat):');
  disp(h);

  zst = z_statistic(P1, P2, D1, D2);
  [h, border_low, border_high] = test_equal_means_z(zst, 0.05);
  disp('(Part 4.3) Error of 2nd kind (Zstat):');
  disp(h);

end

function part5()

  N = 10000;
  D1 = 3;
  D2 = 3;

  sys10 = 0;
  sys20 = 0;

  sys11 = 2;
  sys21 = 0;

  a = 0.05;

  nmin = 3;
  nmax = 40;

  for n = nmin:nmax
    i = n - nmin + 1;
    [hz0(1,i), hz1(1,i), ht0(1,i), ht1(1,i)] = ...
      calculate_all_tests(n, N, D1, D2, sys10, sys20, sys11, sys21, a);
  end

  hz0 = hz0.*(-1) + 1;
  ht0 = ht0.*(-1) + 1;

  figure;
  subplot(2,1,1);
  title('Errors (n)');
  hold on;
  n = nmin:nmax;
  plot(n, hz0, 'DisplayName','Error 1, Z-stats');
  plot(n, ht0, 'DisplayName','Error 1, T-stats');
  plot(n, hz1, 'DisplayName','Error 2, Z-stats');
  plot(n, ht1, 'DisplayName','Error 2, T-stats');
  legend;

  sys11min = 0;
  sys11max = 6;
  n = 8;

  i = 0;
  for sys11 = sys11min:0.1:sys11max
    i = i + 1;
    [hz0(1,i), hz1(1,i), ht0(1,i), ht1(1,i)] = ...
      calculate_all_tests(n, N, D1, D2, sys10, sys20, sys11, sys21, a);
  end

  hz0 = hz0.*(-1) + 1;
  ht0 = ht0.*(-1) + 1;

  subplot(2,1,2);
  title('Errors (delta)');
  hold on;
  sys11 = sys11min:0.1:sys11max;
  plot(sys11, hz0, 'DisplayName','Error 1, Z-stats');
  plot(sys11, ht0, 'DisplayName','Error 1, T-stats');
  plot(sys11, hz1, 'DisplayName','Error 2, Z-stats');
  plot(sys11, ht1, 'DisplayName','Error 2, T-stats');
  legend;

end

function e = e_statistic(P1, P2)
  e = mean(P1-P2);
end

function z = z_statistic(P1, P2, D1, D2)
  n = size(P1,1);
  e = e_statistic(P1, P2);
  z = e./sqrt((D1 + D2)/n);
end

function t = t_statistic(P1, P2)
  e = e_statistic(P1, P2);
  t = e./(std(P1-P2,0)/sqrt(size(P1,1)));
  t = t(find(t > -4 & t < 8));
end

function [h, border_low, border_high] = test_equal_means_e(e, D1, D2, n, alpha)
  borders = norminv([alpha/2, 1-alpha/2], 0, sqrt((D1 + D2)/n));
  border_low = borders(1);
  border_high = borders(2);
  h = length(find(e > border_low &  e < border_high))/length(e);
end

function [h, border_low, border_high] = test_equal_means_z(z, alpha)
  borders = norminv([alpha/2, 1-alpha/2], 0, 1);
  border_low = borders(1);
  border_high = borders(2);
  h = length(find(z > border_low &  z < border_high))/length(z);
end

function [h, border_low, border_high] = test_equal_means_t(t, n, alpha)
  borders = tinv([alpha/2, 1-alpha/2], n-1);
  border_low = borders(1);
  border_high = borders(2);
  h = length(find(t > border_low &  t < border_high))/length(t);
end

% calculates "test equal means"
% of Z and T statistics
% in case of zero and first hypothesis
function [hz0, hz1, ht0, ht1] =  calculate_all_tests...
  (n, N, D1, D2, sys10, sys20, sys11, sys21, a)

  [P10,P20] = generate_sensor_data(n, N, sys10, D1, sys20, D2);
  [P11,P21] = generate_sensor_data(n, N, sys11, D1, sys21, D2);

  z0 = z_statistic(P10, P20, D1, D2);
  z1 = z_statistic(P11, P21, D1, D2);
  t0 = t_statistic(P10, P20);
  t1 = t_statistic(P11, P21);

  [hz0, bl, bh] = test_equal_means_z(z0, a);
  [hz1, bl, bh] = test_equal_means_z(z1, a);
  [ht0, bl, bh] = test_equal_means_t(t0, n, a);
  [ht1, bl, bh] = test_equal_means_t(t1, n, a);

end

% genarates aperiodic pressure data
% with noise and systematic error
function [P1, P2] = generate_sensor_data...
  (n, N, systematic1, dispersion1, systematic2, dispersion2)
  % single aperiodic timeseries
  [t,p] = aperiodic(n);
  % copy timeseries into N columns
  p = repmat(p,1,N);
  % add noise to each column and define P1 and P2
  P1 = p + randn(n,N) * sqrt(dispersion1) + systematic1;
  P2 = p + randn(n,N) * sqrt(dispersion2) + systematic2;
end

% generates aperiodic pressure data
function [t,p] = aperiodic(n)
  t = (1:n)';
  tau = (1/n)*3;
  p = (exp(t.*(-1).*tau).*(-1) + 1).*20 + 10;
end

% histdensity
function [centers, density] = hist_density(x, bin_count)
  dx = (max(x) - min(x)) / bin_count;
  [counts, centers] = hist(x, bin_count);
  density = (counts / length(x)) / dx;
end
