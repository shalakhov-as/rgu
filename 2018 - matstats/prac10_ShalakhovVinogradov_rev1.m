function main()

  clc;
  close all;
  
  % ------ PART 1-4 -------
  
  [u, y] = load_tf_ident_data(4);
  
  n = 2;
  t = (0:length(u)-1)';
  
  H = ident_tf(u, y, n);
  yeval = lsim(H,u-u(1),t);
  
  % ------ PART 5-6 -------
  
  figure;

  subplot(3,1,1);
  hold on;
  title('Input');
  plot(t, u);
  
  subplot(3,1,2);
  hold on;
  title('Responces');
  plot(t, y-y(1));
  plot(t, yeval);
  legend('y data', 'y eval');
  
  subplot(3,1,3);
  hold on;
  title('Responce error');
  plot(t, y-y(1)-yeval);
  
  % ------ PART 7 -------
  
  dydata = diff(y,n);
  dyeval = diff(yeval,n);
  dyerror = dydata - dyeval;
  
  figure;

  subplot(3,1,1);
  hold on;
  title('Input');
  plot(t, u);
  
  subplot(3,1,2);
  hold on;
  title('n-order deriatives');
  plot(dydata);
  plot(dyeval);
  legend('d(n)y data', 'd(n)y eval');
  
  subplot(3,1,3);
  hold on;
  title('Deriatives error');
  plot(dyerror);

end

function tfun = ident_tf(u, y, n)

  y = y - y(1);
  u = u - u(1);
  
  M = build_plan_matrix(u, y, n);
  Y = (diff(y,n));

  b = inv(M' * M) * M' * Y;

  numerator = b(n+1,1);
  denominator = [1, b(1:n,1)'];

  tfun = tf(numerator, denominator);
  
end

function M = build_plan_matrix(u, y, n)
  
  for i = 1:n-1
    deriv = diff(y,i);
    M(:,n-i) = (deriv(1:length(y)-n,1))' .* (-1);
  end
  M(:,n) = (y(1:length(y)-n))' .* (-1);
  M(:,n+1) = (u(1:length(y)-n))';

end

function [u, y] = load_tf_ident_data(var_number)

  fname = ['tf_ident_data/' num2str(var_number) '.csv'];
  data = load(fname);

  u = data(2:end, 1);
  y = data(2:end, 2);

end