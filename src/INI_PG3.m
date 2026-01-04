% written on 04/16/2018

% this is a K-means clustering algorithm for initial values in the mixture
% models. 

% Specifically, the procedure is as follows:

% 1. For input x1, ..., xn, specify the number of clusters K, 
%    place centroids c1, ..., cK at random locations
% 2. For each point xi, 
%        (i) find nearest centroid, i.e., argmin_j  D(xi,
%            cj), where D represents euclidean distance. 
%        (ii) assign xi to cluster j
% 3. For each cluster j = 1, ... K, cj(a) = 1/nj sum_xi belongs to cj
%    xi(a), i.e., new centroid cj = mean of all points xi assigned to
%    cluster j in previous step
% 4. Stop when none of the cluster assignments change. 

% for reference in matlab: https://www.mathworks.com/help/stats/kmeans.html#bues5gz

% for k-means code in matlab, https://www.mathworks.com/help/stats/kmeans.html



% partition the observations of the n-by-p data matrix X into k clusters, and returns an n-by-1 vector (idx) containing cluster indices of each observation. 
% Rows of X correspond to points and columns correspond to variables.

% This is for PG model

% COMMENT: the order is an issue.

function [ini_theta] = INI_PG3(F)

k = 3; % k is the number of clusters
FF = F(:);
idx = kmeans(FF, k); % idx is the index.
% F1 
G = [FF idx];
% F3 are the samples that belong to the first population
F3 = G(G(:, 2)==3);
% F1 are the samples that belong to the second population
F1 = G(G(:, 2)==2);
% F2 are the samples that belong to the third population
F2 = G(G(:, 2)==1);

% For correlated PG model, the mean is alpha/beta, and the variance is
% alpha/beta^2, we may use this to determine the estimate of alpha1, beta1,
% alpha2 and beta2. 

% Another method is to use robust estimate for the mean, i.e., for PG
% model, the mean alpha/beta = median(X), and alpha/beta^2 =
% 1.48*median(X-median(X)),

% We will use the more robust estimator. 

% The proportion estimate is length(F1)/(length(F1)+length(F2)), i.e., 

pi1 = length(F1)/length(FF);
pi2 = length(F2)/length(FF);

% robust estimate of sample variance
% var1 = (1.48*median(F1-median(F1)))^2;
% var2 = (1.48*median(F2-median(F2)))^2;

var1 = var(F1);
var2 = var(F2);
var3 = var(F3);

% robust initial value for PG model
% beta1 = median(F1)/var1;
% alpha1 = median(F1)*beta1;
% 
% beta2 = median(F2)/var2;
% alpha2 = median(F2)*beta2;

beta1 = mean(F1)/var1;
alpha1 = mean(F1)*beta1;

beta2 = mean(F2)/var2;
alpha2 = mean(F2)*beta2;

beta3 = mean(F3)/var3;
alpha3 = mean(F3)*beta3;

% initial value for PL model
% for derivation, see my paper mixture model.
% mu1 = log(median(F1));
% sigma1 = sqrt(log((1+sqrt(1+4*var1/median(F1)^2))/2));
% 
% mu2 = log(median(F2));
% sigma2 = sqrt(log((1+sqrt(1+4*var2/median(F2)^2))/2));

% if alpha1 <= alpha2
%     pi = 1-pi;
%     aa = alpha1;
%     alpha1 = alpha2;
%     alpha2 = aa;
%     bb = beta1;
%     beta1 = beta2;
%     beta2 = bb;
% end

ini_theta = [pi1 alpha1 beta1 pi2 alpha2 beta2 alpha3 beta3];

display(ini_theta)

end






