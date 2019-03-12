function main()
  clc;
  clear all;
  close all;
  warning off;

  % original beta vector
  bi = [272; 0; -1.1*(10^(-5))];
  % dispersion
  Deps = 5;
  n = 20;

  Qmin = 2000;
  Qmax = 3000;

  N = 10000;

  Qplan = part1(n, bi, Qmin, Qmax, Deps);
  bmc = part21(Qplan, bi, Deps, N);

  figure('NumberTitle', 'off', 'Name', 'Distributions of forecast metrics');
  [Qn, Yi, Yevalmc, Dy] = part22(Qplan, bi, bmc, Deps, N);
  part23(Qn, N, bi, Yi, Yevalmc, Dy);
  part24(Qn, N, bi, Yi, Yevalmc, Dy, Deps);

end

% PART 1 - regression
function Qplan = part1(n, bi, Qmin, Qmax, Deps)

  Qplan = build_plan(n, Qmin, Qmax);
  Qn = build_plan(n, 0.8*Qmin, 1.2*Qmax);

  % PART 1

  [Y, X] = generate_process_data(Qplan, bi, Deps);
  beval = ols(Y, X);

  disp('Beta original:');
  disp(bi);
  disp('Beta evaluation:');
  disp(beval);

  [Yi, Xi] = generate_process_data_no_noise(Qn, bi);
  [Yeval, Xeval] = generate_process_data_no_noise(Qn, beval);

  figure('NumberTitle', 'off', 'Name', 'Measurements and characteristics');
  hold on;
  plot(Qn,Yi,'--','LineWidth', 2,'DisplayName','True characteristic');
  plot(Qplan,Y,'o','LineWidth', 2,'DisplayName','Measured datapoints');
  plot(Qn,Yeval,'-','LineWidth', 2,'DisplayName','Evaluated characteristic');
  legend('show');

end

% PART 2.1 - beta distribution
function bmc = part21(Qplan, bi, Deps, N)

  % Monte-Carlo generation
  for i=1:N
    % Ymc is matrix, each column is noised measures
    [Ymc(:,i), X] = generate_process_data(Qplan, bi, Deps);
    % bi is matrix,
    % each column is corresponding beta evaluation
    beval = ols(Ymc(:,i), X);
    bmc(:,i) = beval;
  end

  % calculating theoretical dispersion of beta
  [Y, X] = generate_process_data(Qplan, bi, Deps);
  K = Deps * inv(transpose(X)*X);
  Db = diag(K)

  figure('NumberTitle', 'off', 'Name', 'Distributions of OLS-evaluations');

  % plotting beta distributions
  for i = 1:length(bi)

    subplot(3,1,i);
    hold on;
    legend('show');
    title(['Beta_' num2str(i-1) ' distributions'])

    % beta histdensity
    [b,p] = hist_density (bmc(i,:), 20);
    plot(b, p, '-','LineWidth', 1,'DisplayName', 'histdensity');
    % beta theoretical (normal) distribution
    p = normpdf(b,bi(i), sqrt(Db(i)));
    plot(b, p, '--','LineWidth', 2,'DisplayName','Theoretical');
    % beta confidence interval
    borders = confidence_interval(bi(i), sqrt(Db(i)), 0.05);
    plot(borders, [0 0], 'or', 'DisplayName', 'Confidence interval');
    % beta original value
    plot(bi(i), 0, 'd', ...
      'MarkerSize', 8, 'DisplayName', 'Original value', 'LineWidth', 2);

  end

end

% PART 2.2 - Y forecast distribution
function [Qn, Yi, Yevalmc, Dy] = part22(Qplan, bi, bmc, Deps, N)

  % calculating forecast theoretical dispersions
  % in points Q1 = 0.8 Qmin, Q2 = 1.2 Qmax
  Qn = [0.8*Qplan(1) 1.2*Qplan(end)];
  k = length(bi) - 1;
  X = build_plan_matrix(Qplan, k);
  K = Deps * inv(transpose(X)*X);
  Xn = build_plan_matrix(Qn, length(bi) - 1);
  for i = 1:length(Qn)
    Dy(i) = Xn(i,:) * K * transpose(Xn(i,:))
  end

  % monte-carlo for forecast
  for i = 1:N
    % N predictions
    Yevalmc(:,i) = generate_process_data_no_noise(Qn, bmc(:,i));
  end

  % original values
  Yi = generate_process_data_no_noise(Qn,bi);

  % plotting forecast distributions
  for i = 1:length(Qn)

    subplot(2,3,3*i-2);
    hold on;
    title(['Y forecast dispersion, Q = Q_' num2str(i)]);
    legend('show');

    % forecast histdensity
    [y,p] = hist_density(Yevalmc(i,:), 20);
    plot(y,p,'DisplayName','histdensity');
    % forecast theoretical (normal) distribution
    p = normpdf(y,Yi(i), sqrt(Dy(i)));
    plot(y,p,'--','LineWidth', 2, 'DisplayName', 'theoretical');
    % forecast confidence interval
    borders = confidence_interval(Yi(i), sqrt(Dy(i)), 0.05);
    plot(borders, [0 0], 'or', 'DisplayName', 'Confidence interval');
    % Y original value
    plot(Yi(i), 0, 'd', ...
      'MarkerSize', 8, 'DisplayName', 'Original value', 'LineWidth', 2);

  end

end

% PART 2.3 - True forecast error distributions
function part23(Qn, N, bi, Yi, Yevalmc, Dy)

  % Monte-Carlo for true error
  for i = 1:N
    ei(:,i) = Yevalmc(:,i) - Yi;
  end

  % plotting ei distributions
  for i = 1:length(Qn)

    subplot(2,3,3*i-1);
    hold on;
    title(['Ei dispersion, Q = Q_' num2str(i)]);
    legend('show');

    % ei histdensity
    [y,p] = hist_density(ei(i,:), 20);
    plot(y,p,'DisplayName','histdensity');
    % ei theoretical (normal) distribution
    p = normpdf(y,0, sqrt(Dy(i)));
    plot(y,p,'--','LineWidth', 2, 'DisplayName', 'theoretical');
    % error confidence interval
    borders = confidence_interval(0, sqrt(Dy(i)), 0.05);
    plot(borders, [0 0], 'or', 'DisplayName', 'Confidence interval');
    % e expected value
    plot(0, 0, 'd', ...
      'MarkerSize', 8, 'DisplayName', 'Expected value', 'LineWidth', 2);

  end

end

% PART 2.4 - Forecast measured error distributions
function part24(Qn, N, bi, Yi, Yevalmc, Dy, Deps)

  % Monte-Carlo for measured error
  for i = 1:N
    Yn = generate_process_data(Qn, bi, Deps);
    en(:,i) = Yevalmc(:,i) - Yn;
  end

  % plotting en distributions
  for i = 1:length(Qn)

    subplot(2,3,3*i);
    hold on;
    title(['En dispersion, Q = Q_' num2str(i)]);
    legend('show');

    % en histdensity
    [y,p] = hist_density(en(i,:), 20);
    plot(y,p,'DisplayName','histdensity');
    % en theoretical (normal) distribution
    p = normpdf(y,0, sqrt(Dy(i) + Deps));
    plot(y,p,'--','LineWidth', 2, 'DisplayName', 'theoretical');
    % error confidence interval
    borders = confidence_interval(0, sqrt(Dy(i) + Deps), 0.05);
    plot(borders, [0 0], 'or', 'DisplayName', 'Confidence interval');
    % e expected value
    plot(0, 0, 'd', ...
      'MarkerSize', 8, 'DisplayName', 'Expected value', 'LineWidth', 2);

  end

end

function Qplan = build_plan(n, Qmin, Qmax)
  dQ = (Qmax - Qmin) / (n - 1);
  Qplan = (Qmin:dQ:Qmax)';
end

function X = build_plan_matrix(Qplan, k)
  for i = 0:k
    X(:,i+1) = Qplan.^i;
  end
end

function [Y, X] = generate_process_data(Qplan, bi, Deps)
  [Y,X] = generate_process_data_no_noise(Qplan, bi);
  Y = Y + randn(length(Y),1).* sqrt(Deps);
end

function [Y, X] = generate_process_data_no_noise(Qplan, bi)
  X = build_plan_matrix(Qplan, length(bi)-1);
  for i = 1:length(bi)
    Y(:,i) = X(:,i).*bi(i);
  end
  Y = sum(Y,2);
end

% OLS calculation
function b = ols (Y, X)
  b = inv(transpose(X) * X) * transpose(X) * Y;
end

function [centers, density] = hist_density(x, bin_count)
  dx = (max(x) - min(x)) / bin_count;
  [counts, centers] = hist(x, bin_count);
  density = (counts / length(x)) / dx;
end

function borders = confidence_interval(mu, sigma, alpha)
  borders = norminv([alpha/2, 1-alpha/2], mu, sigma);
end
