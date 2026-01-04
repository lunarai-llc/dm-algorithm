% written on 02/24/2018

% this is to generate iid Y directly from PG model



function [Y] = Generate_Y_iid_PG2(true_theta, n)

% true_pi = 0.3; true_alpha1 = 10; true_beta1 = 1;
% true_alpha2 = 1; true_beta2 = 2;
% 
% n = 1000; 
% 
% true_theta = [true_pi true_alpha1 true_beta1 true_alpha2 true_beta2];

true_pi = true_theta(1); 
true_alpha1 = true_theta(2);
true_beta1 = true_theta(3); 
true_alpha2 = true_theta(4);
true_beta2  = true_theta(5);

Y = zeros(n, 1);

for i = 1:n
    uu = rand;
    if uu <= true_pi
        lambda = gamrnd(true_alpha1, 1/true_beta1);
        Y(i)   = poissrnd(lambda);
    else
        lambda = gamrnd(true_alpha2, 1/true_beta2);
        Y(i)   = poissrnd(lambda);
    end

end


end

