% created by Eisa Keramati
% RBF neural network for classification using ES algorithm
% May 22nd 2019
input = [1 2; 1 3; 4 2];
y_star = [1; 1; 0];
gamma = 1/(10);

%%
clear, clc, close all
fun = 1;                    % function to minimize
switch fun
  case 1
    %  bi-class classification
    %   Minimum: 0
    %f      = @(x,u) (1-x(1,:)).^2 + 100*(x(2,:)-x(1,:).^2).^2;
    %f       = @(x,u) (1/(2*length(input)))*sum(abs((sign(RBF(input, x, gamma, y_star)) - y_star))); 
    f      = @(x,u) (1-x(1,:)).^2 + 100*(x(2,:)-x(1,:).^2).^2;
    %f       = @(x,u) (1/(2))*sum(abs((sign(RBF(input, x, gamma, y_star)) - y_star)));
    dim = 2;
    v_num = 7;
    n_x    = (dim+1)*v_num;                           % 'n_x' states
    limits = repmat([-10 20], n_x, 1);      % Boundaries
    obj    = 0;
otherwise
    error('Not supported equation');
end

%% Setting initial parameters
nf      = 1;                 % length of the output vector 'f(x,y)'
mu      = 100;               % parent population size
lambda  = 100;               % offspring population size
gen     = 10;               % number of generations
sel     = '+';               % Selection scheme (Pag. 78 in (BACK))
rec_obj = 3;                 % Type of recombination to use on object
                             % variables (Pag. 74 in (BACK))
                             % See 'recombination.m'
rec_str = 2;                 % Type of recombination to use on strategy
                             % parameters (Pag. 74 in (BACK))
u        = 0;                % external excitation
func_sel = 3;                % fitness function type selection

%% Run "Evolutionary Strategy" (ES):
[min_x, min_f, off, EPS,idx, arr1, arr2, arr3] = evolution_strategy(f, mu, lambda, gen, sel, rec_obj, rec_str, u, obj, nf, n_x, limits, func_sel, dim);
%disp("start")
%disp(min_x)
%disp("mid")
%disp(min_f)
%disp("done")
disp(arr1)
disp('------------------------------')
disp(arr2)
disp('------------------------------')
disp(arr3)

arr2(arr2 == 0 ) = NaN;
[val,Ind] = min(arr2);
%%    
%v = [1 2.5; 4 2.5];
%res = GMat_calculator(input, v, gamma);
%disp(res);
%res2 = w_calculator(res, y_star);
%disp(res2);
%res3 = y_calculator(res, res2);
%disp(res3);
a = randn(20, 2) + [1 300];
c = randn(20, 2) + [60 300];
b = randn(20, 2) + [400 1];

if (func_sel == 3)
dataTable = readtable('41.csv', 'Format', '%f%f%f');
%dataTable = dataTable(randperm(size(dataTable,1)),:);
dataTable.Properties.VariableNames = {'col1', 'col2', 'col3'};
a = dataTable.col1(1:4000);
b = dataTable.col2(1:4000);
input = [a,b];
y_star = dataTable.col3(1:4000);

end

if (func_sel == 2)
    dataTable = readtable('reg.csv', 'Format', '%f%f%f%f');
%dataTable = dataTable(randperm(size(dataTable,1)),:);
dataTable.Properties.VariableNames = {'col1', 'col2', 'col3', 'col4'};
a = dataTable.col1(1:1500);
%a = normalize(a);
b = dataTable.col2(1:1500);
%b = normalize(b);
c = dataTable.col3(1:1500);
%c = normalize(c);
input = [a,b,c];
y_star = zeros(1500, 1);
for i=1:1:1500
    y_star(i) = dataTable.col4(i);
end
end

if(func_sel == 1)
    dataTable = readtable('22.csv', 'Format', '%f%f%f');
%dataTable = dataTable(randperm(size(dataTable,1)),:);
dataTable.Properties.VariableNames = {'col1', 'col2', 'col3'};
a = dataTable.col1(1:4000);
b = dataTable.col2(1:4000);
input = [a,b];
y_star = dataTable.col3(1:4000);
end

g = GMat_calculator(input, arr1(:,Ind), 1, dim);
if (func_sel==2)
y = g*arr3(:, Ind);
acc = (1/2)*transpose((y - y_star))*(y - y_star);
scatter((1:1:1500),y)
hold('on')
scatter((1:1:1500),y_star)
end

if (func_sel==3)
y = g*arr3(:, :, Ind);
acc = sum(sign((abs(indmax(y) - y_star))))/(length(input));
data = zeros(4000, 2);
data2 = zeros(4000, 2);
a =0;
for i = 1:1:4000
    
    disp(indmax(y(i)))
    disp(y(i))
    disp('----------------------------------------')
    if(indmax(y(i,:)) ~= y_star(i))
        a = a+1;
        disp(indmax(y(i)))
        temp = [input(i,1), input(i,2)];
        data(i,:) = temp;
    else
        temp = [input(i,1), input(i,2)];
        data2(i,:) = temp;
    end
    
end

disp('++++++++++++++++')
v = zeros(v_num, 2);
for i = 1:3:n_x
    temp1 = arr1(:,Ind);
    temp = [temp1(i), temp1(i+1)];
    v(i,:) = temp;
end

scatter(data2(:,1), data2(:,2),'filled')
hold('on')
scatter(data(:,1), data(:,2),'filled')
hold('on')
scatter(v(:,1), v(:,2))

end


if(func_sel == 1)
    dataTable = readtable('22.csv', 'Format', '%f%f%f');
%dataTable = dataTable(randperm(size(dataTable,1)),:);
dataTable.Properties.VariableNames = {'col1', 'col2', 'col3'};
a = dataTable.col1(1:4000);
b = dataTable.col2(1:4000);
input = [a,b];
y_star = dataTable.col3(1:4000);
end

g = GMat_calculator(input, arr1(:,Ind), 1, dim);
if(func_sel == 3)
y = g*arr3(:, :, Ind);
acc = sum(sign((abs(indmax(y) - indmax(y_star)))))/(length(input));
end
if (func_sel==2)
y = g*arr3(:, Ind);
acc = (1/2)*transpose((y - y_star))*(y - y_star);
scatter((1:1:1500),y)
hold('on')
scatter((1:1:1500),y_star)
end

if (func_sel==1)
y = g*arr3(:, Ind);
acc = (1/(2*length(input)))*sum(abs((sign(y) - y_star)));
data = zeros(4000, 2);
data2 = zeros(4000, 2);
a =0;
for i = 1:1:4000
    
    if(sign(y(i)) ~= y_star(i))
        a = a+1;
        temp = [input(i,1), input(i,2)];
        data(i,:) = temp;
    else
        temp = [input(i,1), input(i,2)];
        data2(i,:) = temp;
    end
    
end

disp('++++++++++++++++')
v = zeros(v_num, 2);
for i = 1:3:n_x
    temp1 = arr1(:,Ind);
    temp = [temp1(i), temp1(i+1)];
    v(i,:) = temp;
end

scatter(data2(:,1), data2(:,2),'filled')
hold('on')
scatter(data(:,1), data(:,2),'filled')
hold('on')
scatter(v(:,1), v(:,2))

end

%disp(arr3(:, :, Ind))
%disp(indmax(y))
%disp(indmax(y_star))
%disp('----------------------')
%disp(y_star)
%disp('---------------------')
%disp(y)
disp(acc)
disp(a)


function[G] = GMat_calculator(x, v, gamma, dim)
    %disp(v)
    %disp(v(1:3))
    %disp('--------------------------')
    G = zeros(length(x), length(v)/(dim+1));
    for i=1:1:length(x)
        for j=1:dim+1:length(v)
            G(i,ceil(j/(dim+1))) = exp((-1)*gamma*norm(v(j:j+dim-1), v(j+dim,:), x(i,:)));
            if isnan(exp((-1)*gamma*norm(v(j:j+dim-1), v(j+dim,:), x(i,:))))
            G(i,ceil(j/(dim+1))) = 0;
            end
        end
    end

end

function[res] = norm(v, p, x)
    t = transpose(v);
    res = (x - t)*transpose(x-t)/(p*p);
end


