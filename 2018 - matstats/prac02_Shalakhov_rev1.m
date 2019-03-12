
% ATTENTION
% Function  read_pump_data is copied into this file.
% The path to CSV files can be corrupted due to encoding.

function main()
  clc;
  close all;
  format longG;

  [Q, H, N] = read_pump_data(2, 1);

  ols_approximation(Q,H,'H(Q)');
  ols_approximation(Q,N,'n(Q)');
end

% perform calculations and visualizing
function ols_approximation(Q, H, mode)

  figure('NumberTitle', 'off', 'Name', 'Approximations');
  hold on;

  % plot experimental data
  plot(Q,H,'o','DisplayName','Experimental');

  % calculate beta vector for Lurie model
  [b,br] = approximate_lurie(Q,H,mode);
  % calculate Q-H points of Lurie model
  Q_plot = transpose(0:100:max(Q)*1.5);
  H_lurie = calc_lurie(Q_plot,b,mode);
  H_lurie_r = calc_lurie(Q_plot,br,mode);
  % calculate vector of residuals
  e_lurie = calc_lurie(Q,b,mode) - H;
  % plot H(Q)
  plot(Q_plot,H_lurie,'--','LineWidth', 2,'DisplayName','Lurie');
  plot(Q_plot,H_lurie_r,'DisplayName','Lurie (r)');

  for k = 0:3
    % calculate beta vector for polynom
    [b,br] = approximate_polynom(Q,H,k)
    % calculate Q-H points of polynom
    H_polynom = polynom(Q_plot,b);
    H_polynom_r = polynom(Q_plot,br);
    % calculate vector of residuals
    e_polynom(:,k+1) = polynom(Q,b) - H;
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

  plot (Q,e_lurie,'.');
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

    plot (Q,e_polynom(:,k+1),'.');
    title(['e (polynom k = ' num2str(k) ')']);

    subplot(5,2,(k+2)*2);
    hist(e_polynom(:,k+1));
    xlim([min(min(e_polynom)) max(max(e_polynom))]);
  end

end

% approximation with functions like
% H(Q) = a - bQ^2
% n(Q) = aQ - bQ^2
function [b,br] = approximate_lurie(Q, H, mode)

  if strcmp(mode,'H(Q)')
    k = 0;
  elseif strcmp(mode,'n(Q)')
    k = 1;
  end

  X(:,1) = Q.^k;
  X(:,2) = Q.^2;

  b = ols(H,X);
  br = regress(H,X);

end

% calculate H/n as
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

% approximation with polynom
% of order k
function [b,br] = approximate_polynom(Q, H, k)
  for i = 1:k+1
    X(:,i) = Q.^(i-1);
  end
  b = ols(H,X);
  br = regress(H,X);
end

% calculate H/n as polynom
% of order k
function H = polynom(Q, b)
  H = zeros(length(Q),1);
  for i = 1:length(b)
    H = H + (Q.^(i-1)).*b(i);
  end
end

% OLS calculation
function b = ols (y, X)
  b = inv(transpose(X) * X) * transpose(X) * y;
end

function [Q, H, N] = read_pump_data(station_number, pump_number)
  filename = ['.\na_csvs\' get_pump_string(station_number, pump_number) '.csv'];
  fid = fopen(filename);
  data = textscan(fid, '%s%s%s', 'delimiter', ';');
  fclose(fid);

  % Convert ',' to '.'
  data = cellfun( @(x) str2double(strrep(x, ',', '.')), data, 'uniformoutput', false);
  data = cell2mat(data);

  Q = data(:, 1);
  H = data(:, 2);
  N = data(:, 3);
end

function s = get_pump_string(station_number, pump_number)
  s = [num2str(station_number, '%02d') '_' num2str(pump_number)];
end
