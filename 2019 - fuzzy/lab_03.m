clc;
close all;
warning on;

flow = 68;
water = 57;

% fuzzy control system
fc = readfis('lab_02.fis');

%% task 1.1

% 'flow is comfortable' membership function
mf = fc.input(2).mf(2); 
y1 = evalmf(flow, mf.params, mf.type);

task11 = 1 - y1

%% task 1.2

% (flow is strong) OR (flow is weak)

% flow is weak
mf = fc.input(2).mf(1); 
y1 = evalmf(flow, mf.params, mf.type);
% flow is strong
mf = fc.input(2).mf(3); 
y2 = evalmf(flow, mf.params, mf.type);

task12_minmax = max(y1, y2)
task12_algebr = y1 + y2 - y1*y2

%% task 1.3

% (water is hot) AND (flow is strong)

% water is hot
mf = fc.input(1).mf(3); 
y1 = evalmf(water, mf.params, mf.type);
% flow is strong
mf = fc.input(2).mf(3); 
y2 = evalmf(flow, mf.params, mf.type);

task13_minmax = min(y1, y2)
task13_algebr = y1*y2

%% task 1.4

% (water is cold) OR NOT [(water is warm) AND (flow is comfort)]

% water is cold
mf = fc.input(1).mf(1); 
y1 = evalmf(water, mf.params, mf.type);
% water is warm
mf = fc.input(1).mf(2); 
y2 = evalmf(water, mf.params, mf.type);
% flow is comfort
mf = fc.input(2).mf(2); 
y3 = evalmf(flow, mf.params, mf.type);

task14_minmax = max(y1, (1 -  min(y2, y3)))
not_y23 = 1 - (y2*y3);
task14_algebr = y1 + not_y23 - y1*not_y23

%% task 2

% fuzzy control system
fc = readfis('lab_03.fis');
% input membership functions
Amf = fc.input(1).mf(1); 
Bmf = fc.input(1).mf(2); 
Cmf = fc.input(1).mf(3);
% x axis values
X = 0:0.01:1;
% evalmf of each value
Alh = evalmf(X, Amf.params, Amf.type);
Blh = evalmf(X, Bmf.params, Bmf.type);
Clh = evalmf(X, Cmf.params, Cmf.type);

figure;

%% task 2.1

% A and B

task_21_minmax = min(Alh, Blh);
task_21_algebr = Alh .* Blh;

subplot(3,1,1);
hold on;
legend;
plot(X, Alh, '--', 'LineWidth', 2, 'DisplayName','A');
plot(X, Blh, '--', 'LineWidth', 2, 'DisplayName','B');
plot(X, Clh, '--', 'LineWidth', 2, 'DisplayName','C');
plot(X, task_21_minmax, 'LineWidth', 2, 'DisplayName','A and B (minmax)');
plot(X, task_21_algebr, 'LineWidth', 2, 'DisplayName','A and B (algebr)');

%% task 2.2

% B and (A or C)

task_22_minmax = min(Blh, max(Alh, Clh));
task_22_algebr = Blh .* ( Alh + Clh - Alh.*Clh);

subplot(3,1,2);
hold on;
legend;
plot(X, Alh, '--', 'LineWidth', 2, 'DisplayName','A');
plot(X, Blh, '--', 'LineWidth', 2, 'DisplayName','B');
plot(X, Clh, '--', 'LineWidth', 2, 'DisplayName','C');
plot(X, task_22_minmax, 'LineWidth', 2, 'DisplayName','B and (A or C) (minmax)');
plot(X, task_22_algebr, 'LineWidth', 2, 'DisplayName','B and (A or C) (algebr)');

%% task 2.3

% B or not (A or B or C)

task_23_minmax = max(Blh, 1 - max(Alh, max(Blh,Clh)));
AorB = Alh + Blh - Alh.*Blh;
notAorBorC = 1 - (AorB + Clh - AorB.*Clh);
task_23_algebr = Blh + notAorBorC - Blh.*notAorBorC;

subplot(3,1,3);
hold on;
legend;
plot(X, Alh, '--', 'LineWidth', 2, 'DisplayName','A');
plot(X, Blh, '--', 'LineWidth', 2, 'DisplayName','B');
plot(X, Clh, '--', 'LineWidth', 2, 'DisplayName','C');
plot(X, task_23_minmax, 'LineWidth', 2, 'DisplayName','B or not (A or B or C) (minmax)');
plot(X, task_23_algebr, 'LineWidth', 2, 'DisplayName','B or not (A or B or C) (algebr)');

%% task 3

mf = task_23_minmax;
centroid = sum(X.*mf)/sum(mf);
bisector = defuzz(X,mf,'bisector');
mom = defuzz(X,mf,'mom');
som = defuzz(X,mf,'som');
lom = defuzz(X,mf,'lom');

plot(centroid, 0, 'o', 'LineWidth', 4, 'MarkerSize',10, 'DisplayName','Centroid');
plot(bisector, 0, 'o', 'LineWidth', 4, 'MarkerSize',10, 'DisplayName','Bisector');
plot(mom, 0, 'o', 'LineWidth', 4, 'MarkerSize',10, 'DisplayName','MoM');
plot(som, 0, 'o', 'LineWidth', 4, 'MarkerSize',10, 'DisplayName','SoM');
plot(lom, 0, 'o', 'LineWidth', 4, 'MarkerSize',10, 'DisplayName','LoM');

%%
%
%
%
%
%
%
%
%
%
%
bp = 1;