function [min_x, min_f, off, EPS,j, min_arr, min_arr2, min_w] = evolution_strategy(f, mu, lambda, gen, sel, rec_obj, rec_str, u, obj, nf, n, limits, func_sel, dim, x_0, sigma)
%%  [min_x, min_f, off, EPS,j] = evolution_strategy(f, mu, lambda, gen, sel, rec_obj, rec_str, u, obj, nf, n, limits, x_0, sigma)
%
%   This function tries to optimize the function 'f' using "evolution
%   strategies". The arguments "x_0" and "sigma" are optional.
%   Taken from "Algorithm 4", Section 2.1.6, pag. 81 in (BACK)
%
%   INPUT DATA:
%
%   - f:       Objective function (handle function: f(x,u))
%   - mu:      Parent population size (positive integer number)
%   - lambda:  Offspring population size (positive integer number)
%   - gen:     Number of generations (positive integer number)
%   - sel:     Selection scheme (Pag. 78 in (BACK)):
%                  * ',' = (mu, lambda)-selection scheme
%                  * '+' = (mu + lambda)-selection scheme
%   - rec_obj: Type of recombination to use on objective variables
%              (Pag. 74 in (BACK)):
%                  * 1   = No recombination
%                  * 2   = Discrete recombination
%                  * 3   = Panmictic discrete recombination
%                  * 4   = Intermediate recombination
%                  * 5   = Panmictic intermediate recombination
%                  * 6   = Generalized intermediate recombination
%                  * 7   = Panmictic generalized intermediate recombination
%   - rec_str: Type of recombination to use on strategy parameters
%              (Pag. 74 in (BACK)).
%   - u:       External excitation (if it does not exist, type 0 (zero))
%   - obj:     Vector with the desired results
%   - nf:      Length of the handle function vector (length(f) x 1 vector)
%   - n:       Length of the vector x_0 (positive integer number)
%   - limits:  Matrix with the limits of the variables (nx2 matrix). The
%              first column is the lower boundary, the second column is
%              the upper boundary.
%   - x_0:     Starting point (optional) (nxmu matrix)
%   - sigma:   Cell with covariance matrices (Optional) (1 x mu cell; each
%              cell has to have an nxn symmetric matrix)
%   -func_sel  chosses the type of fitness function(editet, not present in
%              the first version of algorithm).
%
%   OUTPUT DATA:
%
%   - min_x:   Cell with the parent population, and whose last component
%              minimizes the objective function 'f'
%              vector)
%   - min_f:   Cell with the values of the Objective Function 'f'
%              (length(f) x 1 vector)
%   - off:     Cell with the offspring population in each generation
%   - EPS:     Vector with the minimum error of each generation (gen x 1
%              vector)
%   - j:       Number of iterations the algorithm ran (Final number of
%              generations)
%
%   Bibliography:
%
%   - BACK, Thomas. "Evolutionary algorithms in theory and practice".
%     Oxford University Press. New York. 1996.
%
% -------------------------------------------------------
% | Developed by:   Gilberto Alejandro Ortiz Garcia     |
% |                 gialorga@gmail.com                  |
% |                 Universidad Nacional de Colombia    |
% |                 Manizales, Colombia.                |
% -------------------------------------------------------
% -------------------------------------------------------
% | Edited by:      Eisa keramati                       |
% -------------------------------------------------------
%
%   Date: 20-Sep-2011
%% Beggining
%% Initialization:
if (func_sel==3)
dataTable = readtable('42.csv', 'Format', '%f%f%f');
%dataTable = dataTable(randperm(size(dataTable,1)),:);
dataTable.Properties.VariableNames = {'col1', 'col2', 'col3'};
a = dataTable.col1(1:1200);
b = dataTable.col2(1:1200);
input = [a,b];
y_star = zeros(1200, 4);
for i=1:1:1200
      if (dataTable.col3(i) == 1)
        y_star(i,:) = [1,0,0,0];
       elseif (dataTable.col3(i) == 2)
        y_star(i,:) = [0,1,0,0];
        
        elseif (dataTable.col3(i) == 3)
        y_star(i,:) = [0,0,1,0];
        
        elseif (dataTable.col3(i) == 4)
        y_star(i,:) = [0,0,0,1];
        
    end
end
end
if(func_sel==2)
    dataTable = readtable('reg.csv', 'Format', '%f%f%f%f');
dataTable = dataTable(randperm(size(dataTable,1)),:);
dataTable.Properties.VariableNames = {'col1', 'col2', 'col3', 'col4'};
a = dataTable.col1(1:900);
%a = normalize(a);
b = dataTable.col2(1:900);
%b = normalize(b);
c = dataTable.col3(1:900);
%c = normalize(c);
input = [a,b,c];
y_star = zeros(900, 1);
for i=1:1:900
    y_star(i) = dataTable.col4(i);
end
end


if(func_sel == 1)
    dataTable = readtable('21.csv', 'Format', '%f%f%f');
%dataTable = dataTable(randperm(size(dataTable,1)),:);
dataTable.Properties.VariableNames = {'col1', 'col2', 'col3'};
a = dataTable.col1(1:1200);
b = dataTable.col2(1:1200);
input = [a,b];
y_star = dataTable.col3(1:1200);
end

%a = randn(100, 2) + [1 300];
%c = randn(100, 2) + [60 300];
%b = randn(100, 2) + [400 1];
%input = [a;c;b];
 %disp(input)
%input =[a;b];
%disp(input)
%input = [1 300.5; 1 300;1.2 300.6; 400 1; 400.5 1; 405.5 1.2; 2 304; 1.5 323; 420 2];

%y_star = [1; 1; 1; -1; -1; -1; 1; 1; -1];
gamma = 1;
min_arr = zeros(n, gen);
min_arr2 = zeros(gen, 1);
if func_sel == 3
    min_w = zeros(n/(dim+1), length(y_star(1,:)), gen);
else
    min_w = zeros(n/(dim+1), gen);
end
if ((sel ~= ',') && (sel ~= '+'))
  error('not supported selection scheme')
end
if exist('x_0','var')
  if exist('sigma','var')
    [x_0, sigma, alpha] = validate_sizes(mu, n, limits, x_0, sigma);
  else
    [x_0, sigma, alpha] = validate_sizes(mu, n, limits, x_0);
  end
else
  [x_0, sigma, alpha] = validate_sizes(mu, n, limits);
end
e     = 1e-8;                         % maximum permissible error
min_x = cell(1,gen);                  % allocate space in memory for min_x
min_f = cell(1,gen);                  % allocate space in memory for min_f
off   = cell(1,gen);                  % allocate space to store the offspring population
min_x{1} = x_0;                       % first point
value = zeros(nf,mu);                 % allocate space for function evaluation
for i = 1:mu
  %disp(length(x_0(:,i)))%
  %disp(x_0(:,i))%
  if (func_sel == 1)
    value(:,i) = (1/(2*length(input)))*sum(abs((sign(RBF(input, x_0(:,i), gamma, y_star, dim)) - y_star)));
  elseif(func_sel == 2)
    value(:,i) = (1/2)*transpose((RBF(input, x_0(:,i), gamma, y_star, dim) - y_star))*(RBF(input, x_0(:,i), gamma, y_star, dim) - y_star);
  elseif(func_sel == 3)
    value(:,i) = sum(sign((abs(indmax(RBF(input, x_0(:,i) , gamma, y_star, dim)) - indmax(y_star)))))/(length(input));
  
  end
  %disp((1/(2*length(input)))*sum(abs((sign(RBF(input, x_0(:,i), gamma, y_star)) - y_star))));
  %disp('--------------------------')
  %value(:,i) = f(x_0(:,i),u);
end
min_f{1} = value;                     % first approximation
off{1}   = zeros(n,1);
j      = 1;                           % generations counter
jj     = 0;                           % stagnation counter
eps    = abs(obj - value(1,:));       % initial error
EPS    = zeros(gen,1);                % allocate space in memory for minimum error every generation
EPS(1) = min(eps);
t1 = find(eps==min(eps));
%disp(t1)
t1 = t1(1);
min_arr(:,1) = x_0(:,t1);
%disp(t1(1))
min_arr2(1) = min(eps);
if n == 1
  %% Plot function
  figure
  X = linspace(limits(1,1),limits(1,2),100);
  %Y = f(X,[]);
  if (func_sel == 1)
    Y = (1/(2*length(input)))*sum(abs((sign(RBF(input, X, gamma, y_star, dim)) - y_star)));
  elseif(func_sel == 2)
    Y = (1/2)*transpose((RBF(input, X, gamma, y_star, dim) - y_star))*(RBF(input, X, gamma, y_star, dim) - y_star);
  elseif(func_sel == 3)
    Y = sum(sign((abs(indmax(RBF(input, X , gamma, y_star, dim)) - indmax(y_star)))))/(length(input));
  end
  plot(X,Y);
  hold on;
elseif n == 2
  %% Plot function
  figure
  [X,Y] = meshgrid(linspace(limits(1,1),limits(1,2),100),linspace(limits(2,1),limits(2,2),100));
  %Z = reshape(f([X(:) Y(:)]',[]), 100, 100);
  if (func_sel == 1)
    Z = reshape((1/(2*length(input)))*sum(abs((sign(RBF(input,[X(:) Y(:)]' , gamma, y_star, dim)) - y_star))), 100, 100);
  elseif(func_sel == 2)
    Z = reshape((1/2)*transpose((RBF(input, [X(:) Y(:)]' , gamma, y_star, dim) - y_star))*(RBF(input, [X(:) Y(:)]', gamma, y_star, dim) - y_star),100,100);
  elseif(func_sel == 3)
    Z = reshape(sum(sign((abs(indmax(RBF(input, [X(:) Y(:)]' , gamma, y_star, dim)) - indmax(y_star)))))/(length(input)), 100, 100); 
  end
  contour(X,Y,Z,30,'k');      % Contour plot
  grid on
  xlabel('x','FontSize',16);
  ylabel('f(x)','FontSize',16);
  title('Offspring evolution','FontSize',18);
  hold on;
  pcolor(X,Y,Z);              % is really a SURF with its view set to directly above
  shading interp
end
%% Begin ES
%while ((j < gen) && (min(eps) > e))
while ((j < gen))
  %% Print report:  
  if mod(j,5) == 0
    fprintf('\tGeneration j = %4d,  fitness = %g\n',j,min(eps));
  end;
  
  %% Recombine:
  [xr,sigmar,alphar] = recombination(rec_obj,rec_str,n,mu,lambda,min_x{j},sigma,alpha);
  off{j+1}           = xr;            % offspring population
  
  %% Mutation:
  [xm,sigmam,alpham] = mutation(n,lambda,xr,sigmar,alphar,limits);
  
  %% Evaluation:
  phie = zeros(nf,lambda);
  if (func_sel ~= 3)
    phie_0 = zeros(n/(dim+1), lambda);
  else
    phie_0 = zeros(n/(dim+1), length(y_star(1,:)), lambda);
  end
  for i = 1:lambda
    %phie(:,i) = f(xm(:,i),u);
    if (func_sel ~= 3) 
       [temp, phie_0(:,i)] = RBF(input,xm(:,i) , gamma, y_star, dim);
    else
        [temp, phie_0(:, :, i)] = RBF(input,xm(:,i) , gamma, y_star, dim);
    end
    %disp(phie_0(:,i))
    %disp('**********')
    if (func_sel == 1)
        phie(:,i) = (1/(2*length(input)))*sum(abs((sign(temp) - y_star)));
    elseif(func_sel == 2)
        [phie(:,i)] = (1/2)*transpose((temp - y_star))*(temp - y_star);
    elseif(func_sel == 3)
        phie(:,i) = sum(sign((abs(indmax(temp) - indmax(y_star)))))/(length(input));
    end
  end
  epse = abs(obj - phie(1,:));
  %% set the two arrays:
  t1 = find(epse==min(epse));
  t1 = t1(1);
  min_arr(:,j) = xm(:,t1);
  if (func_sel ~= 3) 
       min_w(:,j) = phie_0(:,t1);
  else
       min_w(:, :, j) = phie_0(:,:,t1);
       
       %disp(phie_0(:,:,t1))
       %disp('zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz')
  end
  %disp(t1)
  min_arr2(j) = min(epse);
  
  %% Selection:
  [min_x{j+1}, sigma, alpha] = selection(sel, mu, lambda, epse, eps, xm, sigmam, min_x{j}, sigma, alpham, alpha);
  %% Store better results:
  value = zeros(nf,mu);               % allocate space for function evaluation
  for i = 1:mu
    %value(:,i) = f(min_x{j+1}(:,i),u);
    if (func_sel == 1)
        value(:,i) = (1/(2*length(input)))*sum(abs((sign(RBF(input,min_x{j+1}(:,i) , gamma, y_star, dim)) - y_star)));
    elseif(func_sel == 2)    
        value(:,i) = (1/2)*transpose((RBF(input, min_x{j+1}(:,i), gamma, y_star, dim) - y_star))*(RBF(input,min_x{j+1}(:,i), gamma, y_star, dim) - y_star);
    elseif(func_sel == 3)
        value(:,i) = sum(sign((abs(indmax(RBF(input, min_x{j+1}(:,i) , gamma, y_star, dim)) - indmax(y_star)))))/(length(input));
    end
  end
  min_f{j+1} = value;                 % next approximation
  eps = abs(obj - value(1,:));        % error
  
  EPS(j+1) = min(eps);
  
  %% Stagnation criterion:
  if (EPS(j) == EPS(j+1))
    jj = jj+1;
  else
    jj = 0;
  end
  
  %% Increase generation counter:
  
  j = j+1;
  
  %% Plot preliminary results
  if n == 1
    %% Plot offspring population at each generation
    h = plot(off{j-1}(1,:),min_f{j-1}(1,:),'*r');
    pause(0.2)
    if (((EPS(j) > e) && (j < gen)) && (jj<30))
      delete(h);
    end
  elseif n == 2
    %% Plot offspring population at each generation
    h = plot(off{j-1}(1,:),off{j-1}(2,:),'*r');
    axis([limits(1,1) limits(1,2) limits(2,1) limits(2,2)]);
    pause(0.2)
    if (((EPS(j) > e) && (j < gen)) && (jj<30))
      delete(h);
    end
  end
  if (jj == 30)
    fprintf('\n\n\tError remains constant for 30 consecutive generations\n\n');
    break
  end
end
end
%% END