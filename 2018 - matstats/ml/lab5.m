function lab5()

  clc;
  clear all;
  close all;


  % PART 1
  a = 4;
  b = 5;
  x1 = 1:100;
  y1 = hygernd(1000,40,50,1,100).*a + b;
  EY = mean(y1)
  SKOY = std(y1)
  DY = SKOY^2

  figure;
  hold on;
  plot(x1,y1,'ro');

  a = 7;
  b = 8;
  x2 = 1:100;
  y2 = hygernd(1000,40,50,1,100).*a + b;
  EY = mean(y2)
  SKOY = std(y2)
  DY = SKOY^2

  plot(x2,y2,'bo');

  figure;
  subplot(1,2,1);

  histogram(y1);
  title('Первое распределение');
  subplot(1,2,2);

  histogram(y2);
  title('Второе распределение');

  %  Гидрослюда в угле  № проб  Зольность, %  Гидрослюда в угле
  data = [1  28.1  7.05  28  10.1  1.96;
  2  15.4  2.62  29  14.0  2.07;
  3  12.8  3.34  30  19.4  6.35;
  4  11.1  1.96  31  13.9  2.42;
  5  10.4  2.37  32  8.1  4.19;
  6  9.3  2.22  33  10.1  1.02;
  7  17.8  2.88  34  10.2  3.73;
  8  13.7  2.53  35  11.8  3.21;
  9  17.8  6.94  36  23.8  4.51;
  10  15.2  3.78  37  9.3  2.70;
  11  16.4  2.61  38  24.5  4.90;
  12  12.8  1.41  39  15.7  3.94;
  13  12.3  1.82  40  16.1  4.92;
  14  13.9  1.78  41  7.5  2.12;
  15  16.1  2.68  42  12.4  2.53;
  16  9.6  3.27  43  13.7  6.43;
  17  6.6  0.57  44  9.0  2.00;
  18  17.6  3.50  45  10.9  4.23;
  19  10.9  1.72  46  21.3  2.34;
  20  8.0  2.20  47  8.3  2.03;
  21  12.1  2.88  48  11.8  3.02;
  22  8.1  2.05  49  24.0  10.52;
  23  13.7  1.60  50  6.9  0.48;
  24  17.6  9.77  51  7.0  3.44;
  25  12.5  3.68  52  47.9  27.67;
  26  15.4  4.55  53  42.3  24.92;
  27  12.5  3.4  54  41.5  9.90;];

  zolnost = [data(:,2); data(:,5)];
  hydroslyuda = [data(:,3); data(:,6)];

  figure;
  plot(zolnost, hydroslyuda,'ro');

  R = corrcoef(zolnost, hydroslyuda)

  fun = @(a, x) a(1) + x.*a(2) + x.^2.*a(3);
  a0 = [0 0 0];
  ac = lsqcurvefit(fun, a0, zolnost, hydroslyuda);

  hold on;
  x = 0:50;
  y = fun(ac, x);
  plot(x,y,'-');

  data2 = [
  232 119 137 25;
  293 120 177 56;
  87 130 98 101;
  121 417 25 634;
  422 355 115 340;
  1580 198 360 195;
  835 567 195 158;
  204 504 493 24;
  218 574 487 210;
  243 404 379 50;
  146 502 247 228;
  49 697 116 335;
  174 579 629 153;
  0 451 260 77;
  0 627 68 195;
  0 597 254 234;
  0 726 211 219;
  0 686 254 75;
  0 683 82 43;
  0 525 100 0;
  0 605 9 0;
  0 1042 30 0;
  0 504 0 0 ;
  0 648 0 0;
  0 220 0 0;]

  rt525 = data2(:,1);
  rt525 = rt525(find(rt525 ~= 0))
  rt518 = data2(:,2);
  rt518 = rt518(find(rt518 ~= 0))
  rt501 = data2(:,3);
  rt501 = rt501(find(rt501 ~= 0))
  rt509 = data2(:,4);
  rt509 = rt509(find(rt509 ~= 0))

  figure;
  boxplot(data2,'notch','on','labels',{'telo525','telo518','telo501','telo509'});

  meanrt525 = mean(rt525)
  meanrt518 = mean(rt518)
  meanrt501 = mean(rt501)
  meanrt509 = mean(rt509)

  freq525 = tabulate(rt525)
  freq518 = tabulate(rt518)
  freq501 = tabulate(rt501)
  freq509 = tabulate(rt509)

  figure;
  subplot(4,1,1);
  histogram(rt525);
  subplot(4,1,2);
  histogram(rt518);
  subplot(4,1,3);
  histogram(rt501);
  subplot(4,1,4);
  histogram(rt509);

  SKO525 = std(rt525)
  SKO518 = std(rt518)
  SKO501 = std(rt501)
  SKO509 = std(rt509)

end
