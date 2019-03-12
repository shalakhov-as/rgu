function prac11()

  clc;
  clear all;
  close all;
  
  n = 15;
  bi = [265; -0.045];
  bapr = [255; -0.055];
  
  R = diag(rand(n,1));
  Q = diag(rand(length(bi),1));
  
  [x,y] = gendata(bi, n, 10);
  
  mmnk(x,y,bi,bapr,n,R,Q);
  ols_with_lims(x,y,bi,bapr,n,R,Q);

end

% PART 1 - MMNK
function bmmnk = mmnk(x,y,bi,bapr,n,R,Q)

  X = build_plan_matrix(x,length(bi)-1);
  
  K = inv(Q + X'*R*X)*X'*R;
  bmmnk = bapr + K*(y - X*bapr);
  disp('MMNK:');
  disp(bmmnk);
  
  [b0,b1,S] = mmnk_quality_data(X,y,R,Q,bapr);
  Smin = mmnk_quality(bmmnk,X,y,R,Q,bapr);
  
  figure;
  hold on;
  title('MMNK minimisation function I(beta)')
  su = surf(b0,b1,log(S));
  su.EdgeColor = 'none';
  scatter3(bmmnk(1),bmmnk(2),log(Smin),'ro','MarkerFaceColor','r');
  legend('I(beta)', 'I(mmnk)');
  
  figure;
  hold on;
  title('MMNK minimisation function I(beta)')
  contour(b0,b1,log(S));
  scatter3(bmmnk(1),bmmnk(2),log(Smin),'ro','MarkerFaceColor','r');
  legend('I(beta)', 'I(mmnk)');
  
end

% PART 2 - MNK(OLS) with limitations
function bqlim = ols_with_lims(x,y,bi,bapr,n,R,Q)

  X = build_plan_matrix(x,length(bi)-1);
  
  H = (X'*R*X + Q)*2;
  f = -2*(y'*R*X + bapr'*Q);
  
  bmin = [-inf; -0.03];
  bmax = [240; +inf];
  
  bquad = quadprog(H,f);
  bqlim = quadprog(H, f, [], [], [], [], bmin, bmax);
  disp('QUADPROG:');
  disp(bquad);
  disp('QUADPROG W/LIMITS:');
  disp(bqlim);
  
  [b0,b1,S] = mmnk_quality_data(X,y,R,Q,bapr);
  Squad = mmnk_quality(bquad,X,y,R,Q,bapr);
  Sqlim = mmnk_quality(bqlim,X,y,R,Q,bapr);
  
  figure;
  hold on;
  title('Quadprog evaluation results')
  contour(b0,b1,log(S));
  plot(bquad(1),bquad(2),'ro','MarkerFaceColor','r');
  plot(bqlim(1),bqlim(2),'bo','MarkerFaceColor','b');
  plot([150 400],[bmin(2) bmin(2)]);
  plot([bmax(1) bmax(1)],[-0.07 -0.01]);
  legend('I(beta)',...
         'Quadprog beta',...
         'Quadprog w/limits',...
         'Limit 1',...
         'Limit 2');
  
end

% generates experimantal data with noise
function [x,y] = gendata(bi, n, D)

  x = 2000:10:(2000 + 10*n - 1);
  bm = fliplr(bi');
  
  y = polyval(bm, x);
  y = y + randn(1,length(x)).*sqrt(D);
  
  x = x';
  y = y';
  
end

% calculates mmnk quality with given beta
function I = mmnk_quality(b,X,Y,R,Q,bapr)
  I = (Y - X*b)'*R*(Y - X*b) + (b - bapr)'*Q*(b - bapr);
end

% calculates mmnk quality with generated
% array of beta values (for surfaces)
function [b0,b1,S] = mmnk_quality_data(X,Y,R,Q,bapr)

  % defining meshgrid around "true" values
  b0_values = 150:1:400;
  b1_values = -0.07:0.0001:-0.01;
  [b0, b1] = meshgrid(b0_values, b1_values);
  
  % calculating minimisation function
  for i = 1:size(b0,1)
    for j = 1:size(b0,2)
      b_current = [b0(i,j); b1(i,j)];
      S(i,j) = mmnk_quality(b_current,X,Y,R,Q,bapr);
    end
  end
  
end

% builds plan matrix
function X = build_plan_matrix(x, k)
  for i = 0:k
    X(:,i+1) = x.^i;
  end
end