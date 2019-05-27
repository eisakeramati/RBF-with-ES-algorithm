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
    n_x    = 4;                           % 'n_x' states
    limits = repmat([0 410], n_x, 1);      % Boundaries
    obj    = 0;
otherwise
    error('Not supported equation');
end

%% Setting initial parameters
nf      = 1;                 % length of the output vector 'f(x,y)'
mu      = 100;               % parent population size
lambda  = 100;               % offspring population size
gen     = 100;               % number of generations
sel     = '+';               % Selection scheme (Pag. 78 in (BACK))
rec_obj = 2;                 % Type of recombination to use on object
                             % variables (Pag. 74 in (BACK))
                             % See 'recombination.m'
rec_str = 2;                 % Type of recombination to use on strategy
                             % parameters (Pag. 74 in (BACK))
u       = 0;                 % external excitation
func_sel = 2;                % fitness function type selection

%% Run "Evolutionary Strategy" (ES):
[min_x, min_f, off, EPS,idx, arr1, arr2] = evolution_strategy(f, mu, lambda, gen, sel, rec_obj, rec_str, u, obj, nf, n_x, limits, func_sel);
%disp("start")
%disp(min_x)
%disp("mid")
%disp(min_f)
%disp("done")
disp(arr1)
disp('------------------------------')
disp(arr2)

%%    
%v = [1 2.5; 4 2.5];
%res = GMat_calculator(input, v, gamma);
%disp(res);
%res2 = w_calculator(res, y_star);
%disp(res2);
%res3 = y_calculator(res, res2);
%disp(res3);


