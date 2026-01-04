% written on 04/12/2018

% this is to generate iid 3 component mixture Y directly from PG model 

function [Y] = Generate_Y_iid_PG3(true_theta, n)

% true_pi1 = 0.5; true_alpha1 = 10; true_beta1 = 1; true_pi2 = 0.25;
% true_alpha2 = 5; true_beta2 = 2; true_alpha3 = 1; true_beta3 = 2;
% 
% n = 10000; 
% 
% true_theta = [true_pi1 true_alpha1 true_beta1 true_pi2 true_alpha2 true_beta2 true_alpha3 true_beta3];

true_pi1 = true_theta(1); 
true_alpha1 = true_theta(2);
true_beta1 = true_theta(3);

true_pi2 = true_theta(4);
true_alpha2 = true_theta(5);
true_beta2  = true_theta(6);

true_alpha3 = true_theta(7);
true_beta3 = true_theta(8);

Y = zeros(n, 1);

for i = 1:n
    uu = rand;
    if uu <= true_pi1
        lambda = gamrnd(true_alpha1, 1/true_beta1);
        Y(i)   = poissrnd(lambda);
    elseif true_pi1 < uu <= (true_pi1 + true_pi2)
        lambda = gamrnd(true_alpha2, 1/true_beta2);
        Y(i)   = poissrnd(lambda);
    else 
        lambda = gamrnd(true_alpha3, 1/true_beta3);
        Y(i)   = poissrnd(lambda);
    end
end





end

