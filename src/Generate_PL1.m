% written on 02/24/2018

% this is to generate iid Y directly from PL model

function [Y] = Generate_PL1(true_theta, n)

% true_mu = 2; true_sigma = .5;
% 
% n = 200;
% 
% true_theta = [true_mu true_sigma];

true_mu = true_theta(1);
true_sigma = true_theta(2); 

Y = zeros(n, 1);

for i = 1:n    
    lambda = lognrnd(true_mu, true_sigma);        
    Y(i)   = poissrnd(lambda);
end


end

