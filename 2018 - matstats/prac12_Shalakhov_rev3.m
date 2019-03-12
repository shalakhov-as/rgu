function prac12()
  clc;
  close all;
  
  [Qe, He, Ne] = read_pump_data(1,1);
  
  n = length(Qe);
  
  b = approximate_lurie(Qe,He);
  Hm = calc_lurie(Qe,b);
  
  e2m = (He - Hm).^2;

  
  Deps = sum(e2m)/n;
  
  Se = sum(e2m)/(n-2);

  e2l = (loo(Qe,He,Se)).^2;
  Dloo = sum(e2l)/n;
  
  disp('D(LOO):');
  disp(Dloo);
  disp('e^2 mean value:')
  disp(Deps);
  disp('S2e:')
  disp(Se);

end

function e = loo(Qe,He,Se)

  hits = 0;
  e = zeros(length(Qe),1);

  for i = 1:length(Qe)
    
    Ql = Qe;
    Hl = He;
    % leaving i-element out
    Ql(i) = [];
    Hl(i) = [];
    % calculatung beta evaluations
    bl = approximate_lurie(Ql,Hl);
    
    % calculating standart deviation
    Hm = calc_lurie(Ql,bl);
    % calculating examination value
    Hexam = calc_lurie(Qe(i),bl);
    
    a = 0.05;
    X = build_plan_matrix(Qe);
    De = Se * (X(i,:)* inv(X'*X)*X(i,:)'+1);
    dtHexam = (Hexam - He(i)) / sqrt(De);
    tbounds = tinv([a/2, 1 - a/2],length(Qe)-length(bl));
    
    % check if experimantal value
    % fits into confidence bounds
    if (dtHexam >= tbounds(1) & dtHexam <= tbounds(2))
       hits = hits + 1;
    end
    
    % save approximation error in array
    e(i) = Hexam - He(i);
    
  end

  percent = hits/length(Qe)*100;
  disp('Hits in 95% interval (%)');
  disp(percent);
    
end

% approximation with
% H(Q) = a - bQ^2
function [b] = approximate_lurie(Q, H)

  X(:,1) = zeros(length(Q),1) + 1;
  X(:,2) = Q.^2;

  b = regress(H,X);

end

function X = build_plan_matrix(Q)

  X(:,1) = zeros(length(Q),1) + 1;
  X(:,2) = Q.^2;

end

% calculate H
% H(Q) = a - bQ^2
function H = calc_lurie(Q, b)
  H = (Q.^0).*b(1) + (Q.^2).*b(2);
end

function [Q, H, N] = read_pump_data(station_number, pump_number)
  filename = ['.\na_csvs\' get_pump_string(station_number, pump_number) '.csv'];
  fid = fopen(filename);
  data = textscan(fid, '%s%s%s', 'delimiter', ';');
  fclose(fid);

  % Convert ',' to '.'
  data = cellfun( @(x) str2double(strrep(x, ',', '.')), data, 'uniformoutput', false);
  data = cell2mat(data);

  Q = data(:, 1);
  H = data(:, 2);
  N = data(:, 3);
end

function s = get_pump_string(station_number, pump_number)
  s = [num2str(station_number, '%02d') '_' num2str(pump_number)];
end
