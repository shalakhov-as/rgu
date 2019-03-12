function prac13()
  clc;
  clear all;
  close all;
  
  filename = 'models_Shalakhov.xlsx';
  e = [ 1; 2; 3; -4; -5; -6; 7; 8; 9; -9; 9; -9; -9; 8];
  n = get_series_count(e);
  [Q, H, N] = read_pump_data(2, 1);
  [beta_qh, S2h, R2h, R2hregress, Rh_analysis] = ols_approximation(Q, H, 'H(Q)');
  [beta_qn, S2n, R2n, R2nregress, Rn_analysis] = ols_approximation(Q, H, 'n(Q)');
  
  column_names = ["b0", "b1", "b2", "b3", ...
    "R2", "R2 (regress)", "Residuals analysis", "Dispersion evalualion"];
  xlswrite(filename,column_names,1,'B2');
  xlswrite(filename,beta_qh,1,'B3');
  xlswrite(filename,R2h,1,'F3');
  xlswrite(filename,R2hregress,1,'G3');
  xlswrite(filename,Rh_analysis,1,'H3');
  xlswrite(filename,S2h,1,'I3');
end



%% Dispersion calculation
function S2 = calc_residual_dispersion(y, y_prediction, k)
  S2 = sum((y - y_prediction).^2) / (length(y) - k);
end

%% Total SS calculation
function TSS = calc_TSS(y)
  TSS = sum((y - mean(y)).^2);
end

%% Residuals SS calculation
function RSS = calc_RSS(y, y_prediction)
  RSS = sum((y - y_prediction).^2);
end

%% Explained SS calculation
function ESS = calc_ESS(y, y_prediction)
  TSS = calc_TSS(y);
  RSS = calc_RSS(y, y_prediction);
  ESS = TSS - RSS;
end

%% R2 calculation
function R2 = calc_R2(y, y_prediction)
  TSS = calc_TSS(y);
  RSS = calc_RSS(y, y_prediction);
  R2 = 1 - (RSS/TSS);
end

%% Sort vector of residuals in order of x
function e_sorted = sort_e_following_x(e,x)
  [xs, ind] = sort(x);
  e_sorted = e(ind);
end

%% Calculationg amount of series
function Ns = get_series_count(e)
  e_shifted = [-e(1); e(1:length(e)-1)];
  e_products = e .* e_shifted;
  [row,col,v] = find(e_products < 0);
  Ns = length(v);
end

%% Count signed elements
function [Nplus, Nminus] = get_signed_values_count(e)
  [r,c,v] = find(e >= 0);
  Nplus = length(v);
  Nminus = length(e) - Nplus;
end

%% Calculate randomness stats
function [random_residuals, z, Ns] = test_residuals_randomness(e)
  Ns = get_series_count(e);
  [Np, Nm] = get_signed_values_count(e);
  
  ENs = 2*Np*Nm / (Np + Nm) + 1;
  DNs = (2*Np*Nm / ((Np + Nm)^2)) * ((2*Np*Nm - (Np + Nm)) / (Np + Nm - 1));
  
  z = (Ns - ENs)/sqrt(DNs);
  pz = normpdf(z,0,1);
  a = 0.05;
  
  if (pz > a/2)
    random_residuals = 0;
  else
    random_residuals = 1;
  end
  
  [h,p,stats] = runstest(e);
  disp('Z-stats of e-randomness test difference:');
  disp(z);
  disp('Z-stats of e-randomness test difference (runstest):');
  disp(stats.z);
end

%% perform OLS calculations and visualizing
function [beta_table, S2, R2, R2regress, R_analysis] = ols_approximation(Q, H, mode)

  figure('NumberTitle', 'off', 'Name', 'Approximations');
  hold on;

  % plot experimental data
  plot(Q,H,'o','DisplayName','Experimental');

  % calculate beta vector for Lurie model
  [b,br,R2l] = approximate_lurie(Q,H,mode);
  beta_table(1,:) = [b(1) 0 b(2) 0];
  R2regress(1,:) = R2l;
  disp('LURIE');
  disp(br);
  % calculate Q-H points of Lurie model
  Q_plot = transpose(0:100:max(Q)*1.5);
  H_lurie = calc_lurie(Q_plot,b,mode);
  H_lurie_r = calc_lurie(Q_plot,br,mode);
  % calculate vector of residuals
  Hl = calc_lurie(Q,b,mode);
  e_lurie = Hl - H;
  e_lurie_ordered = sort_e_following_x(e_lurie,Q);
  R_analysis(1,:) = test_residuals_randomness(e_lurie_ordered);
  S2(1,:) = calc_residual_dispersion(H, Hl, 2);
  R2(1,:) = calc_R2(H, Hl);
  % plot H(Q)
  plot(Q_plot,H_lurie,'--','LineWidth', 2,'DisplayName','Lurie');
  plot(Q_plot,H_lurie_r,'DisplayName','Lurie (r)');

  for k = 0:3
    % calculate beta vector for polynom
    [b,br,R2p] = approximate_polynom(Q,H,k)
    beta_table(k+2,1:k+1) = b;
    R2regress(k+2,:) = R2p;
    % calculate Q-H points of polynom
    H_polynom = polynom(Q_plot,b);
    H_polynom_r = polynom(Q_plot,br);
    % calculate vector of residuals
    Hp = polynom(Q,b);
    e_polynom(:,k+1) = Hp - H;
    e_polynom_ordered = sort_e_following_x(e_polynom(:,k+1),Q);
    R_analysis(k+2,:) = test_residuals_randomness(e_polynom_ordered);
    S2(k+2,:) = calc_residual_dispersion(H, Hp, k);
    R2(k+2,:) = calc_R2(H, Hp);
    % plot H(Q)
    plot(Q_plot,H_polynom,'--','LineWidth', 2,'DisplayName',['Polynom k = ', num2str(k)]);
    plot(Q_plot,H_polynom_r,'DisplayName',['Polynom k = ', num2str(k), ' (r)']);
  end

  title(mode);
  legend('show');
  legend('boxoff');
  xlim([0 max(Q)*1.5]);
  ylim([0 max(H)*1.5]);

  % visualize residuals data
  % e(Q) and histograms
  figure('NumberTitle', 'off', 'Name', [mode ' residuals']);
  subplot(5,2,1);
  hold on;

  ylim([min(min(e_polynom)) max(max(e_polynom))]);
  xlim([min(Q) max(Q)]);
  line(xlim, [0 0], 'Color','black');  %x-axis
  
  [Q,ind] = sort(Q);
  y = e_lurie(ind);
  plot (Q,y,'.');
  title('e (Lurie)');

  subplot(5,2,2);
  hist(e_lurie);
  xlim([min(min(e_polynom)) max(max(e_polynom))]);

  for k = 0:3
    subplot(5,2,(k+2)*2-1);
    hold on;

    ylim([min(min(e_polynom)) max(max(e_polynom))]);
    xlim([min(Q) max(Q)]);

    line(xlim, [0 0], 'Color','black');  %x-axis
    
    y = e_polynom(:,k+1);
    y = y(ind);
    plot (Q,y,'.');
    title(['e (polynom k = ' num2str(k) ')']);

    subplot(5,2,(k+2)*2);
    hist(e_polynom(:,k+1));
    xlim([min(min(e_polynom)) max(max(e_polynom))]);
  end

end

%% approximation with functions like
% H(Q) = a - bQ^2
% n(Q) = aQ - bQ^2
function [b,br,R2] = approximate_lurie(Q, H, mode)

  if strcmp(mode,'H(Q)')
    k = 0;
  elseif strcmp(mode,'n(Q)')
    k = 1;
  end

  X(:,1) = Q.^k;
  X(:,2) = Q.^2;

  b = ols(H,X);
  [br,ccc,ddd,fff,stats] = regress(H,X);
  R2 = stats(1);

end

%% calculate H/n as
% H(Q) = a - bQ^2
% n(Q) = aQ - bQ^2
function H = calc_lurie(Q, b, mode)

  if strcmp(mode,'H(Q)')
    k = 0;
  elseif strcmp(mode,'n(Q)')
    k = 1;
  end

  H = (Q.^k).*b(1) + (Q.^2).*b(2);
end

%% approximation with polynom
% of order k
function [b,br,R2] = approximate_polynom(Q, H, k)
  for i = 1:k+1
    X(:,i) = Q.^(i-1);
  end
  b = ols(H,X);
  [br,ccc,ddd,fff,stats] = regress(H,X);
  R2 = stats(1);
end

%% calculate H/n as polynom
% of order k
function H = polynom(Q, b)
  H = zeros(length(Q),1);
  for i = 1:length(b)
    H = H + (Q.^(i-1)).*b(i);
  end
end

%% OLS calculation
function b = ols (y, X)
  b = inv(transpose(X) * X) * transpose(X) * y;
end