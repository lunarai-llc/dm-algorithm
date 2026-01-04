% this is the example on how to use this Discrete kernel package

%% Discrete Model Example
% Usage: est_theta = Discrete_Disparity(X, disparity, model); 

% X: input data

% disparity: disparity method: 'hd'(Hellinger distance method), 
% 'ned' (Negative Exponential Disparity Method)

% model: 'binom' (Binomial model), 'poiss' (Poisson model), 
% 'geo' (Geometric model), 'nb' (Negative Binomial model), 
% 'pg' (Poisson-Gamma model), 'pl'(Poisson-Lognormal model) 

% generate 100 Poisson random variables with mean 3
X = poissrnd(3, 100, 1);
% plot histogram
histogram(X);
% Choice #1: use Discrete_Disparity funtion hd method to estimate the parameter
hd_theta = Discrete_Disparity(X, 'hd', 'poiss');
display(hd_theta)

% Choice #2: use Discrete_Disparity funtion ned method to estimate the parameter
ned_theta = Discrete_Disparity(X, 'ned', 'poiss');
display(ned_theta)



