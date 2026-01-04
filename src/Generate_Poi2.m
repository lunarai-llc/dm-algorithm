% written on 04/20/2018


% this is to generate poisson mixture.

% pi = 0.2;
% 
% n = 500;
% 
% lambda = [1 2];
% 
% true_theta = [pi lambda];

function [Y] = Generate_Poi2(true_theta, n)

pi = true_theta(1);

lambda = zeros(1, 2);
lambda(1) = true_theta(2);
lambda(2) = true_theta(3);

% generate mixture poisson 
Y = zeros(n, 1);
for i = 1:n
    uu = rand;
    if (uu <= pi)
        Y(i) = poissrnd(lambda(1));
    else
        Y(i) = poissrnd(lambda(2));
    end        
end


end
