% written 04/21/2018

% this is initial value for 2-component poisson mixture


function [ini_theta] = INI_Poi2(F)

k = 2; % k is the number of clusters
FF = F(:);
idx = kmeans(FF, k); % idx is the index.
% F1 
G = [FF idx];
% F1 are the samples that belong to the first population
F1 = G(G(:, 2)==2);
% F2 are the samples that belong to the second population
F2 = G(G(:, 2)==1);

% For correlated PG model, the mean is alpha/beta, and the variance is
% alpha/beta^2, we may use this to determine the estimate of alpha1, beta1,
% alpha2 and beta2. 

% Another method is to use robust estimate for the mean, i.e., for PG
% model, the mean alpha/beta = median(X), and alpha/beta^2 =
% 1.48*median(X-median(X)),

% We will use the more robust estimator. 

% The proportion estimate is length(F1)/(length(F1)+length(F2)), i.e., 

pi = length(F1)/length(FF);

lambda1 = mean(F1);
lambda2 = mean(F2);

if lambda1 >= lambda2
    aa = lambda1;
    lambda1 = lambda2;
    lambda2 = aa;
    pi = 1-pi;
end

% if alpha1 <= alpha2
%     pi = 1-pi;
%     aa = alpha1;
%     alpha1 = alpha2;
%     alpha2 = aa;
%     bb = beta1;
%     beta1 = beta2;
%     beta2 = bb;
% end

ini_theta = [pi lambda1 lambda2];

display(ini_theta)

end

