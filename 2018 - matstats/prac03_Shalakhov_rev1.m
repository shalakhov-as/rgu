
% ATTENTION
% Function  read_pump_data is copied into this file.
% The path to CSV files can be corrupted due to encoding.

function main()
  clc;
  close all;

  % "true" values
  b = [336.228736243678; -0.0534243863519426];
  % generating data sample including noise
  [Q, H] = generate_qh_exp_data(b, 3, 15);
  % visualize data sample
  figure;
  plot(Q,H,'o');

  % defining meshgrid around "true" values
  b0_values = 300:5:400;
  b1_values = (-0.06:0.001:-0.05);
  [b0, b1] = meshgrid(b0_values, b1_values);

  % calculating sum of differences' squares
  % for each beta combination
  for i = 1:size(b0,1)
    for j = 1:size(b0,2)
      b_current = [b0(i,j) b1(i,j)];
      S(i,j) = log(sum_of_squares(Q,H,b_current));
    end
  end

  [b_eval br_eval] = approximate_polynom(Q, H, 1)
  s_eval = log(sum_of_squares(Q,H,b_eval));
  sr_eval = log(sum_of_squares(Q,H,br_eval));

  % visualize surface of diffs' sums
  figure;
  hold on;
  surf(b0,b1,S);
  scatter3(b_eval(1),b_eval(2),s_eval, 'o','MarkerFaceColor','r');
  scatter3(br_eval(1),br_eval(2),sr_eval, 'o','MarkerFaceColor','g');

  % visualize levels
  figure;
  hold on;
  contour(b0,b1,S);
  scatter3(b_eval(1),b_eval(2),s_eval, 'o','MarkerFaceColor','r');
  scatter3(br_eval(1),br_eval(2),sr_eval, 'o','MarkerFaceColor','g');
end

% generator of dispersed data
function [Q, H] = generate_qh_exp_data(b, sko, count)
  Q = (2000:(1000/count):3000)';
  H = (Q.^0).*b(1) + (Q.^1).*b(2) + randn(length(Q),1)*sko;
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

% OLS calculation
function b = ols (y, X)
  b = inv(transpose(X) * X) * transpose(X) * y;
end

% calculating target function
function S = sum_of_squares(Q,H,b)
  H_eval = QH_linear(Q, b);
  E = H_eval - H;
  S = sum(E.^2);
end

function H = QH_linear(Q, b)
  H = (Q.^0).*b(1) + (Q.^1).*b(2);
end
