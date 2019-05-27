
function[y] = RBF(x,v, gamma, y_star)
    y_0 = GMat_calculator(x,v,gamma);
    %disp(y_0)
    %disp('---------------------')
    y_1 = w_calculator(y_0, y_star);
    y = y_calculator(y_0, y_1);
    %disp(y)
    %disp('---------------------')
end

function[G] = GMat_calculator(x, v, gamma)
    G = zeros(length(x), length(v)/2);
    for i=1:1:length(x)
        for j=1:2:length(v)
            G(i,ceil(j/2)) = exp((-1)*gamma*norm(v(j,:), (v(j+1,:)), x(i,:)));
            if isnan(exp((-1)*gamma*norm(v(j,:), (v(j+1,:)), x(i,:))))
            G(i,ceil(j/2)) = 0;
            end
        end
    end

end

function[res] = norm(v, w, x)
    t = [transpose(v) transpose(w)];
    res = (x - t)*transpose(x-t);
end

function[w] = w_calculator(G, y)
    G_trans = transpose(G);
    %eigenvalue = eig(G_trans*G, 'matrix');
    inverse = inv(G_trans*G + 5*eye(length(G_trans*G)));
    w = inverse*G_trans*y;
end

function[y] = y_calculator(G, W)
    y = G*W;
end



%function[error] = Regression_error(y, y_star)
%    error = (1/2)*transpose((y-y_star))*(y-y_star);
%end

%function[error] = Biclassification_error(y, y_star)
%    error = 0;
%    for i = 1 : 1 : length(y)
%        error = error + abs((sign(y(i)) - y_star(i)));
%    end
%    error = 1 - (1/(2*length(y)))*error;
%end
