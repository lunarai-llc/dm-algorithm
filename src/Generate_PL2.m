% written on 04/20/2018

% this is to generate iid Y directly from PL model

% this is exactly the same as Generate_Y_iid_PL(true_theta, n)

function [Y] = Generate_PL2(true_theta, n)

% true_pi = 0.7; true_mu1 = 2; true_sigma1 = .5;
% true_mu2 = 1; true_sigma2 = 1;
% 
% n = 200;
% 
% true_theta = [true_pi true_mu1 true_sigma1 true_mu2 true_sigma2];

true_pi = true_theta(1); true_mu1 = true_theta(2);
true_sigma1 = true_theta(3); 
true_mu2 = true_theta(4);
true_sigma2  = true_theta(5);

% lambda1 = lognrnd(true_mu1, true_sigma1, n, 1);
% 
% lambda2 = lognrnd(true_mu2, true_sigma2, n, 1);
% 
% Y = zeros(n, 1);
% 
% for i = 1:n    
%     uu = rand;
%     if uu <= true_pi
%         Y(i) = poissrnd(lambda1(i));
%     else
%         Y(i) = poissrnd(lambda2(i));
%     end    
% end


Y = zeros(n, 1);

% for i = 1:n
%     uu = rand;
%     if uu <= true_pi
%         lambda = lognrnd(true_mu1, true_sigma1);
%         Y(i)   = poissrnd(lambda);
%     else
%         lambda = lognrnd(true_mu2, true_sigma2);
%         Y(i)   = poissrnd(lambda);
%     end
% 
% end


for i = 1:n
    uu = rand;
    if uu <= true_pi
        lambda = lognrnd(true_mu1, true_sigma1);
        % Y(i)   = poissrnd(lambda);
    else
        lambda = lognrnd(true_mu2, true_sigma2);
        % Y(i)   = poissrnd(lambda);
    end
    Y(i)   = poissrnd(lambda);
end


end

