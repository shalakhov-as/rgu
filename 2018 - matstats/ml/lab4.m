function lab4()

  clc;
  close all;
  clear all;

  % PART 1

  % a
  F = @(x) 1./sqrt(25 + 3 * x);
  I = quad(F, 0, -3)

  % b
  F = @(x) sinh(x).*sinh(x);
  I = quad(F, 0, 1)

  % c
  F = @(x,y) x * sin(y) + y * sin(x);
  I = dblquad(F, 1, 2, 0, 1)

  % PART 2

  S = [30 40 50 60 70 80 90 100 110 120];
  G = [6.36 6.85 7.34 7.84 8.08 8.32 8.57 8.7 8.82 8.94];

  fun = @(x,S) x(1) + x(2).*S + x(3).*S.^2;
  x0 = [0 0 0]

  x = lsqcurvefit(fun,  x0, S, G)

  Slin = linspace(S(1), S(end));
  figure;
  plot(S, G, 'ro', Slin, fun(x,Slin));

  eps = sqrt(sum((G - fun(x, S)).^2)./length(S))
  ymin = fun(x, S(1));
  delta = (eps / ymin) * 100

  % PART 3

  P3S = dsolve('D2y + 6*Dy - 11*y = sin(t)')

  [t,y] = ode45(@part3eq, [0 20], [1; 0]);

  figure;
  plot(t, y(:,1), '--o', t, y(:,2), '--o')

  xlabel('time');
  ylabel('y');
  legend('y1','y2');

  ylim([0 10]);

  % PART 4

  % a
  syms x
  D1 = diff((2*x+3)/(x*x-5*x+5))
  D2 = diff(D1)
  % b
  syms x
  D1 = diff(1/x + 2*log(x) - log(x)/x)
  D2 = diff(D1)

  % PART 5

  % a
  syms x
  f = 1 / (x^4 * sqrt(4 - x^2));
  I5a = int(f)

  % b
  syms x
  f = (1 + sqrt(x)) / x^2;
  I5b = int(f, [1 4])

  % c
  syms x
  f = x * exp(x);
  I5c = int(f, [-inf 0])

  % d
  syms x
  f = 1 / ((x^2 + 1)^2);
  I5d = int(f, [-inf +inf])

  % e
  syms x
  f = atan(x) / ((x^2 + 1));
  I5e = int(f, [0 +inf])

  % f
  syms x
  f = x / ((1 - x^2));
  I5f = int(int(f))

  % PART 6

  % a
  s6a = dsolve('D2y - 2*Dy = x^2 - 1')
  % b
  s6b = dsolve('D2y - 2*Dy + 10*y = sin(3*x) + exp(x)')

end

function dydt = part3eq(t, y)
  dydt = [y(2); -6*y(2) + 11*y(1) + sin(t)];
end
