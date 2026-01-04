% written on 04/12/2018

% this is to generate iid 4 component mixture Y directly from PG model 

function [Y] = Generate_Y_iid_PG5(true_theta, n)

% true_pi1 = 0.30;  true_alpha1 = 20; true_beta1 = 1.5; 
% true_pi2 = 0.25;  true_alpha2 = 10; true_beta2 = 1; 
% true_pi3 = 0.20;  true_alpha3 = 5;  true_beta3 = 2;
% true_pi4 = 0.15;  true_alpha4 = 3;  true_beta4 = 1.5;
%                   true_alpha5 = 1;  true_beta5 = 2;
% 
% n = 10000; 
% 
% true_theta = [true_pi1 true_alpha1 true_beta1  ...
%               true_pi2 true_alpha2 true_beta2  ...
%               true_pi3 true_alpha3 true_beta3  ...
%               true_pi4 true_alpha4 true_beta4  ...
%                        true_alpha5 true_beta5  ];
                
true_pi1 = true_theta(1); 
true_alpha1 = true_theta(2);
true_beta1 = true_theta(3);

true_pi2 = true_theta(4);
true_alpha2 = true_theta(5);
true_beta2  = true_theta(6);

true_pi3 = true_theta(7);
true_alpha3 = true_theta(8);
true_beta3 = true_theta(9);

true_pi4 = true_theta(10);
true_alpha4 = true_theta(11);
true_beta4 = true_theta(12);

true_alpha5 = true_theta(13);
true_beta5 = true_theta(1);

Y = zeros(n, 1);

for i = 1:n
    uu = rand;
    if uu <= true_pi1
        lambda = gamrnd(true_alpha1, 1/true_beta1);
        Y(i)   = poissrnd(lambda);
    elseif true_pi1 < uu <= (true_pi1 + true_pi2)
        lambda = gamrnd(true_alpha2, 1/true_beta2);
        Y(i)   = poissrnd(lambda);
    elseif (true_pi1 + true_pi2)  < uu <= (true_pi1 + true_pi2 + true_pi3)
        lambda = gamrnd(true_alpha3, 1/true_beta3);
        Y(i)   = poissrnd(lambda);
    elseif (true_pi1 + true_pi2 + true_pi3)  < uu <= (true_pi1 + true_pi2 + true_pi3 + true_pi4)
        lambda = gamrnd(true_alpha4, 1/true_beta4);
        Y(i)   = poissrnd(lambda);    
    else
        lambda = gamrnd(true_alpha5, 1/true_beta5);
        Y(i)   = poissrnd(lambda);    
    end
end



end

