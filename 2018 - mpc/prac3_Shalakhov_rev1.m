function prac3_Shalakhov_rev1()

  clc;
  clear;
  close all;
  warning off all;
  
  disp('Revision 1 covers parts 1,2, and 3.1');
  disp('MPC with applied constraints will be covered in revision 2');
  disp('___');

  % timeline definitions
  dt = 0.05;
  stop_time = 15;
  time = 0:dt:stop_time-dt;

  % transfer function definitions
  gain = 10;
  roots = [-2, -1.2 + 2i, -1.2 - 2i];
  p = poly(roots);
  numerator = gain * p(length(p));
  transfer_function = tf(numerator, p);
  [A, B, C] = space_state_model_disc(transfer_function, dt);
  
  disp('calculating part 1 ...');
  prac3_part1(A,B,C,time,gain);
  disp('part 1 done, calculating part 2 ...');
  prac3_part2(A,B,C,time);
  disp('part 2 done, calculating part 3.1 ...');
  prac3_part3_1(A,B,C,time);
  disp('part 3.1 done');

end

% part 1 - calculation incremental model response 
% and comparison vs model from prac2 
function prac3_part1(A,B,C,time,gain)
  
  [Ai, Bi, Ci] = get_increment_model(A,B,C);

  P = 200;
  L = get_L_matrix(A,C,P);
  Li = get_L_matrix(Ai,Ci,P);
  M = get_M_matrix(A,B,C,P);
  Mi = get_M_matrix(Ai,Bi,Ci,P);
  
  x0 = [0;0;0];
  u = ones(P,1);
  
  p0 = [x0;0];
  v = eye(P,1);
  
  y = calculate_responce(L, M, x0, u);
  z = calculate_responce(Li, Mi, p0, v);
  
  t = time(1:length(y));
  
  figure('Name','Part 1. Models comparison');
  title('Original discrete model vs incremental model');
  hold on;
  plot (t, y, 'Color', [69/255, 104/255, 160/255], 'LineWidth', 2);
  plot (t, z, '.', 'Color', [226/255, 129/255, 160/255], 'MarkerSize', 12);
  plot (t, u, 'Color', [55/255, 183/255, 55/255], 'LineWidth', 2);
  plot (t, u.*gain, '--', 'Color', [229/255, 13/2557, 0/255], 'LineWidth', 2);
  grid on;
  legend('simple discrete model: y(t)', 'incremental discrete model: z(t)', 'u(t)', 'gain level');

end

% part 2 - astatic MPC controller
function prac3_part2(A,B,C,time)

  rk = 2;
  SP = build_sp_matrix(rk,length(time));
  
  X0 = [0;0;0];
  Rk = [1];
  Qk = [1];
  
  % changing P
  P = 8;
  [Y08,U08] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  P = 16;
  [Y16,U16] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  P = 32;
  [Y32,U32] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  
  % changing R,Q
  P = 16;
  Rk = 1;
  Qk = 10;
  [Y1_10,U1_10] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  Rk = 10;
  Qk = 1;
  [Y10_1,U10_1] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  Rk = 1;
  Qk = 1;
  [Y1_1,U1_1] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk);
  
  figure('Name','Part 2. Astatic MPC modelling');
  
  subplot(3,2,1);
  title('MPC control, P = 8');
  hold on;
  grid on;
  plot (time, Y08, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U08, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
  
  subplot(3,2,3);
  title('MPC control, P = 16');
  hold on;
  grid on;
  plot (time, Y16, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U16, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
  
  subplot(3,2,5);
  title('MPC control, P = 32');
  hold on;
  grid on;
  plot (time, Y32, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U32, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
  
  
  subplot(3,2,2);
  title('MPC control, R = 1, Q = 10');
  hold on;
  grid on;
  plot (time, Y1_10, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U1_10, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
  
  subplot(3,2,4);
  title('MPC control, R = 10, Q = 1');
  hold on;
  grid on;
  plot (time, Y10_1, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U10_1, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
  
  subplot(3,2,6);
  title('MPC control, R = 1, Q = 1');
  hold on;
  grid on;
  plot (time, Y1_1, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U1_1, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
end

% part 3.1 - astatic MPC with quadprog without constraints
% comparison with part 2
function prac3_part3_1(A,B,C,time)

  rk = 2;
  SP = build_sp_matrix(rk,length(time));
  
  X0 = [0;0;0];
  Rk = [1];
  Qk = [1];
  
  % changing P
  P = 8;
  [Y08,U08] = mqp_modelling(A,B,C,time,P,rk,X0,Rk,Qk,[],[]);
  P = 16;
  [Y16,U16] = mqp_modelling(A,B,C,time,P,rk,X0,Rk,Qk,[],[]);
  P = 32;
  [Y32,U32] = mqp_modelling(A,B,C,time,P,rk,X0,Rk,Qk,[],[]);
  
  % changing R,Q
  P = 16;
  Rk = 1;
  Qk = 10;
  [Y1_10,U1_10] = mqp_modelling(A,B,C,time,P,rk,X0,Rk,Qk,[],[]);
  Rk = 10;
  Qk = 1;
  [Y10_1,U10_1] = mqp_modelling(A,B,C,time,P,rk,X0,Rk,Qk,[],[]);
  Rk = 1;
  Qk = 1;
  [Y1_1,U1_1] = mqp_modelling(A,B,C,time,P,rk,X0,Rk,Qk,[],[]);
  
  figure('Name','Part 3. Quadprog MPC modelling');
  
  subplot(3,2,1);
  title('MPC control, P = 8');
  hold on;
  grid on;
  plot (time, Y08, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U08, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
  
  subplot(3,2,3);
  title('MPC control, P = 16');
  hold on;
  grid on;
  plot (time, Y16, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U16, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
  
  subplot(3,2,5);
  title('MPC control, P = 32');
  hold on;
  grid on;
  plot (time, Y32, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U32, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
  
  
  subplot(3,2,2);
  title('MPC control, R = 1, Q = 10');
  hold on;
  grid on;
  plot (time, Y1_10, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U1_10, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
  
  subplot(3,2,4);
  title('MPC control, R = 10, Q = 1');
  hold on;
  grid on;
  plot (time, Y10_1, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U10_1, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');
  
  subplot(3,2,6);
  title('MPC control, R = 1, Q = 1');
  hold on;
  grid on;
  plot (time, Y1_1, 'LineWidth', 2, 'Color', [69/255, 104/255, 160/255]);
  plot (time, U1_1, 'LineWidth', 2, 'Color', [229/255, 13/2557, 0/255]);
  plot (time, SP, '--', 'LineWidth', 2, 'Color', [226/255, 129/255, 160/255]);
  legend('y(t)', 'u(t)', 'set point');

end

% part 3.2 - astatic MPC with constraints
function prac3_part3_2(A,B,C,time)

  

end

% gets matrixes of incremental model
function [Ai,Bi,Ci] = get_increment_model(A,B,C)

  n = size(A,1);
  r = size(C,1);
  
  Ai = [A, zeros(n,r); C*A, eye(r,r)];
  Bi = [B; C*B];
  Ci = [zeros(r,n), eye(r,r)];

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

% calculates responce of SSM
function y = calculate_responce(L, M, x0, u)
 
  y = L*x0 + M*u;

end

function [Y, U] = mpc_modelling(A,B,C,time,P,rk,X0,Rk,Qk)

  X = [X0,X0];
  Y = [C*X0,C*X0];
  U = [0];
  
  [Ai, Bi, Ci] = get_increment_model(A,B,C);
  Li = get_L_matrix(Ai,Ci,P);
  Mi = get_M_matrix(Ai,Bi,Ci,P);
  
  for k = 2:length(time)
    
    pk = [X(:,k) - X(:,k-1); Y(:,k)];
    
    dU = mpc_controller(pk, P, Li, Mi, Rk, Qk, rk);
    
    U = [U, U(k-1) + dU(1)];
    Xnext = A * X(:,k) + B * U(:, k);
    Ynext = C * X(:,k);
    
    X = [X, Xnext];
    Y = [Y, Ynext];
    
  end
  
  Y = Y(1:end-1);
  
end

% mqp stands for Mpc QuadProg
% Ac, Bc - constraint equation system coefficients
function [Y, U] = mqp_modelling(A,B,C,time,P,rk,X0,Rk,Qk,Ac,Bc)

  X = [X0,X0];
  Y = [C*X0,C*X0];
  U = [0];
  
  [Ai, Bi, Ci] = get_increment_model(A,B,C);
  Li = get_L_matrix(Ai,Ci,P);
  Mi = get_M_matrix(Ai,Bi,Ci,P);
  
  for k = 2:length(time)
    
    pk = [X(:,k) - X(:,k-1); Y(:,k)];
    
    dU = mpc_quadprog(pk, P, Li, Mi, Rk, Qk, rk, Ac, Bc);
    
    U = [U, U(k-1) + dU(1)];
    Xnext = A * X(:,k) + B * U(:, k);
    Ynext = C * X(:,k);
    
    X = [X, Xnext];
    Y = [Y, Ynext];
    
  end
  
  Y = Y(1:end-1);
  
end

function dUk = mpc_controller(Xk, P, L, M, R, Q, rk)
  
  Rk = R;
  Qk = Q;
  R = get_R_matrix(Rk, P);
  Q = get_Q_matrix(Qk, P);
  
  r = build_sp_matrix(rk,P);
  
  T = -inv(Q + (M')*R*M)*(M')*R;
  dUk = T*(L*Xk - r);

end

function dUk = mpc_quadprog(pk, P, L, M, R, Q, rk, Ac, Bc)
  
  Rk = R;
  Qk = Q;
  R = get_R_matrix(Rk, P);
  Q = get_Q_matrix(Qk, P);
  
  r = build_sp_matrix(rk,P);
  
  H = (M')*R*M + Q;
  f = (M')*R*L*pk - (M')*R*r;
  %g = (pk')*(L')*R*L*pk + (r')*R*r - 2*(pk')*(L')*R*r;
  
  options = optimset('Display', 'off');
  dUk = quadprog(H,f,Ac,Bc,[],[],[],[],[],options);
  
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

function R = get_R_matrix(Rk, P)

  R = Rk;
  for i = 2:P
    R = blkdiag(R,Rk);
  end
  
end

function Q = get_Q_matrix(Qk, P)

  Q = get_R_matrix(Qk,P);
  
end

% builds SP matrix [rP x 1] using constant rk
function r = build_sp_matrix(rk,P)

  r = [];
  for i = 1:P
    r = [r; rk];
  end

end