% written on 04/20/2018

% this is to generate iid 3 component mixture Y directly from PL model 

% true_pi1 = 0.15; true_mu1 = 2; true_sigma1 = .5; 
% true_pi2 = 0.2;  true_mu2 = 1; true_sigma2 = .5; 
% true_pi3 = 0.3;  true_mu3 = 3; true_sigma3 = .5;
%                  true_mu4 = 2.5; true_sigma4 = 0.5;
% n = 500; 
% 
% true_theta = [true_pi1 true_mu1 true_sigma1 ...
%               true_pi2 true_mu2 true_sigma2 ...
%               true_pi3 true_mu3 true_sigma3 ...
%                        true_mu4 true_sigma4];
% 
% Y = Generate_PL4(true_theta, n);

function [Y] = Generate_PL4(true_theta, n)

% true_pi1 = 0.15; true_mu1 = 2; true_sigma1 = .5; 
% true_pi2 = 0.2;  true_mu2 = 1; true_sigma2 = .5; 
% true_pi3 = 0.3;  true_mu3 = 3; true_sigma3 = .5;
%                  true_mu4 = 2.5; true_sigma4 = 0.5;
% n = 500; 
% 
% true_theta = [true_pi1 true_mu1 true_sigma1 ...
%               true_pi2 true_mu2 true_sigma2 ...
%               true_pi3 true_mu3 true_sigma3 ...
%                        true_mu4 true_sigma4];

true_pi1 = true_theta(1); 
true_mu1 = true_theta(2);
true_sigma1 = true_theta(3);

true_pi2 = true_theta(4);
true_mu2 = true_theta(5);
true_sigma2  = true_theta(6);

true_pi3 = true_theta(7);
true_mu3 = true_theta(8);
true_sigma3 = true_theta(9);

true_mu4 = true_theta(10);
true_sigma4 = true_theta(11);


Y = zeros(n, 1);

for i = 1:n
    uu = rand;
    if uu <= true_pi1
        lambda = lognrnd(true_mu1, true_sigma1);
        Y(i)   = poissrnd(lambda);
    elseif true_pi1 < uu <= (true_pi1 + true_pi2)
        lambda = lognrnd(true_mu2, true_sigma2);
        Y(i)   = poissrnd(lambda);
    elseif (true_pi1 + true_pi2)  < uu <= (true_pi1 + true_pi2 + true_pi3)
        lambda = lognrnd(true_mu3, true_sigma3);
        Y(i)   = poissrnd(lambda);
    else
        lambda = lognrnd(true_mu4, true_sigma4);
        Y(i)   = poissrnd(lambda);
    end
end

end

