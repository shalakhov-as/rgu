clc;

figure;
subplot(3,1,1);
fis = readfis('lab_01_01.fis');
plotmf(fis, 'input', 1);
grid on;

subplot(3,1,2);
fis = readfis('lab_01_02.fis');
plotmf(fis, 'input', 1);
grid on;

subplot(3,1,3);
fis = readfis('lab_01_03.fis');
plotmf(fis, 'input', 1);
grid on;

figure;

subplot(2,1,1);
fis = readfis('lab_01_04.fis');
plotmf(fis, 'input', 1);
grid on;

subplot(2,1,2);
fis = readfis('lab_01_05.fis');
plotmf(fis, 'input', 1);
grid on;
