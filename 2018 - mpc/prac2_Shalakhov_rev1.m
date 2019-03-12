function prac2_Shalakhov_rev1()

  clc;
  clear;
  close all;

  % timeline definitions
  dt = 0.05;
  stop_time = 10;
  time = 0:dt:stop_time-dt;
  P = length(time);

  % transfer function definitions
  gain = 10;
  roots = [-2, -0.6 + 2i, -0.6 - 2i];
  p = poly(roots);
  numerator = gain * p(length(p));
  transfer_function = tf(numerator, p);
  [A, B, C] = space_state_model_disc(transfer_function, dt);
  
  % executing prac parts
  prac2_part1(A,B,C,time);
  prac2_part2(A,B,C,time);
  
end

function prac2_part1(A,B,C,time)
  
  P = length(time);
  x0 = [0;0;0];
  u = ones(P,1);
  u(1) = 0;
  y = calculate_responce(A, B, C, P, x0, u);
  
  figure;
  title('y = L*x0 * M*u');
  hold on;
  plot (time, y, 'LineWidth', 2);
  plot (time, u, 'LineWidth', 2);
  legend('y(t)', 'u(t)');

end


function prac2_part2(A,B,C,time)

  rk = 3;
  X0 = [0;0;0];
  Rk = [10];
  Qk = [1];
  
  % changing P
  P = 3;
  [Y10,U10] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  P = 7;
  [Y30,U30] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  P = 10;
  [Y90,U90] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  
  % changing R,Q
  P = 6;
  Rk = 1;
  Qk = 10;
  [Y1_10,U1_10] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  Rk = 10;
  Qk = 1;
  [Y10_1,U10_1] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  Rk = 1;
  Qk = 1;
  [Y1_1,U1_1] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  
  figure;
  
  subplot(3,2,1);
  title('MPC control, P = 3');
  hold on;
  grid on;
  plot (time, Y10, 'LineWidth', 2);
  plot (time, U10, 'LineWidth', 2);
  legend('y(t)', 'u(t)');
  
  subplot(3,2,3);
  title('MPC control, P = 7');
  hold on;
  grid on;
  plot (time, Y30, 'LineWidth', 2);
  plot (time, U30, 'LineWidth', 2);
  legend('y(t)', 'u(t)');
  
  subplot(3,2,5);
  title('MPC control, P = 10');
  hold on;
  grid on;
  plot (time, Y90, 'LineWidth', 2);
  plot (time, U90, 'LineWidth', 2);
  legend('y(t)', 'u(t)');
  
  
  subplot(3,2,2);
  title('MPC control, R = 1, Q = 10');
  hold on;
  grid on;
  plot (time, Y1_10, 'LineWidth', 2);
  plot (time, U1_10, 'LineWidth', 2);
  legend('y(t)', 'u(t)');
  
  subplot(3,2,4);
  title('MPC control, R = 10, Q = 1');
  hold on;
  grid on;
  plot (time, Y10_1, 'LineWidth', 2);
  plot (time, U10_1, 'LineWidth', 2);
  legend('y(t)', 'u(t)');
  
  subplot(3,2,6);
  title('MPC control, R = 1, Q = 1');
  hold on;
  grid on;
  plot (time, Y1_1, 'LineWidth', 2);
  plot (time, U1_1, 'LineWidth', 2);
  legend('y(t)', 'u(t)');
end


function [Y, U] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk)

  X = X0;
  Y = C * X0;
  U = [];
  
  for k = 1:length(time)
    
    Uoptk = mpc_controller(X(:,k), P, A, B, C, Rk, Qk, rk);
    U = [U, Uoptk(1)];
    
    Xnext = A * X(:,k) + B * U(:, k);
    Ynext = C * X(:,k);
    
    X = [X, Xnext];
    Y = [Y, Ynext];
    
  end
  
  Y = Y(1:end-1);
  
end

function R = get_R_matrix(Rk, P)

  R = Rk;
  for i = 2:P
    R = blkdiag(R,Rk);
  end
  
end

function Q = get_Q_matrix(Qk, P)

  Q = get_R_matrix(Qk,P);
  
end

function L = get_L_matrix(A,C,P)

  L = C*A;
  for i = 2:P
    L = [L; C*(A^i)];
  end

end

function M = get_M_matrix(A,B,C,P)

  M = [];
  for i = 1:P
    
    line = [];
    
    for j = 1:i
      line = [line, C*(A^(i-j))*B];
    end
    
    for j = (i+1):P
      line = [line, zeros(size(line,1))];
    end
    
    M = [M; line];
  end

end

% builds SP matrix [rP x 1] with constant rk
function r = build_sp_matrix(rk,P)

  r = zeros(0);
  for i = 1:P
    r = [r; rk];
  end

end

% calculates responce of SSM
function y = calculate_responce(A, B, C, P, x0, u)

  L = get_L_matrix(A,C,P);
  M = get_M_matrix(A,B,C,P);
 
  y = L*x0 + M*u;

end

function Uk = mpc_controller(Xk, P, A, B, C, R, Q, rk)
  
  Rk = R;
  Qk = Q;
  R = get_R_matrix(Rk, P);
  Q = get_Q_matrix(Qk, P);
  L = get_L_matrix(A,C,P);
  M = get_M_matrix(A,B,C,P);
  r = build_sp_matrix(rk,P);
  
  T = -inv(Q + (M')*R*M)*(M')*R;
  Uk = T*(L*Xk - r);

end

% space-state model for discrete time (from prac1)
function [A, B, C] = space_state_model_disc(tf, dt)
k = tf.Numerator{1}(4);
p = tf.Denominator{1};

A(1,:) = [-p(2)/p(1)*dt + 1, -p(3)/p(1)*dt, -p(4)/p(1)*dt];
A(2,:) = [dt, 1, 0];
A(3,:) = [0, dt, 1];
B = [k*dt; 0; 0];
C = [0, 0, 1];
end