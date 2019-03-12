function prac1_Shalakhov_rev1()

clc;
close all;
clear all;

% timeline definitions
dt = 0.01;
time = 0:dt:10;

% transfer function definitions
gain = 10;
roots = [-2, -0.6 + 2i, -0.6 - 2i];
p = poly(roots);
numerator = gain * p(length(p));
transfer_function = tf(numerator, p);

% continuous time models calculations
ss_model_matlab = ss(transfer_function);
[Ac, Bc, Cc] = space_state_model_cont(transfer_function);

% discrete time models calculations
ss_model_matlab_disc = c2d(ss_model_matlab, dt);
[Ad, Bd, Cd] = space_state_model_disc(transfer_function, dt);

% part 1 - models comparison

y_tf = step(transfer_function, time);
y_c2d = step(ss_model_matlab_disc, time);
[~, y_manual] = calculate_response(Ad, Bd, Cd, time, 'disable_PID');

figure('Name','Space-state models');
hold on;

title('Space-state models');
plot (time, y_tf, 'LineWidth', 2);
plot (time, y_c2d, ':r', 'LineWidth', 2);
plot (time, y_manual, '--', 'LineWidth', 2);
grid on;
legend('ML continuous time model', 'ML discrete time model', 'Custom discrete time model');

% part 2 - PID control

[~, y, U] = calculate_response(Ad, Bd, Cd, time, 'enable_PID');
U = U(1:length(y));

figure('Name','PID control');
subplot(2,1,1);
hold on;
title('Process variable');
plot (time, y, 'LineWidth', 2);
legend('y(t)');
grid on;
subplot(2,1,2);
title('Control variable');
hold on;
plot (time, U, ':r', 'LineWidth', 2);
legend('U(t)');
grid on;


end

% space-state model for discrete time
function [A, B, C] = space_state_model_disc(tf, dt)
k = tf.Numerator{1}(4);
p = tf.Denominator{1};

A(1,:) = [-p(2)/p(1)*dt + 1, -p(3)/p(1)*dt, -p(4)/p(1)*dt];
A(2,:) = [dt, 1, 0];
A(3,:) = [0, dt, 1];
B = [k*dt; 0; 0];
C = [0, 0, 1];
end

% space-state model for continuous time
function [A, B, C] = space_state_model_cont(tf)
k = tf.Numerator{1}(4);
p = tf.Denominator{1};

A(1,:) = [-p(2)/p(1), -p(3)/p(1), -p(4)/p(1)];
A(2,:) = [1, 0, 0];
A(3,:) = [0, 1, 0];
B = [k; 0; 0];
C = [0, 0, 1];

end

% calculates respose of space-state model
% for each point of time
function [X, Y, U] = calculate_response(A, B, C, time, enable_PID)
X(:,1) = [0;0;0];
U(:,1) = 0;
Kc = 0.07;
SP = 4;
for k = 1:length(time)
    X(:,k+1) = A * X(:,k) + B * U(:,k);
    Y(:,k) = C * X(:,k);
    
    if (strcmp(enable_PID,'enable_PID'))
        U(:,k+1) = pid_controller(Kc, SP - Y(:,k));
    else
        U(:,k+1) = 1;
    end
end
end

function u = pid_controller(Kc, e)
u = Kc*e + 0.4;
end