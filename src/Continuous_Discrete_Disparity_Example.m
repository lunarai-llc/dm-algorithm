% this is the example on how to use this package

%% 1. Discrete Model Example
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
% use Discrete_Disparity funtion hd method to estimate the parameter
hd_theta = Discrete_Disparity(X, 'hd', 'poiss');
display(hd_theta)

% use Discrete_Disparity funtion ned method to estimate the parameter
ned_theta = Discrete_Disparity(X, 'ned', 'poiss');
display(ned_theta)

%% 2. Continuous Model Example
% Usage: Continuous_Disparity(X, model)

% X: input data

% model: continuous model: The choices are: 
% 'normal', 'gamma', 'beta', 'lognormal', 'weibull', 'logistic' ,'t', 
% 'exp', 'rayleigh', 'nakagami', 'rician', 'chisq', 'inversegaussian'

% default disparity method: 'hd'; another choice: 'ned'
% default bandwidth selection: see bandwidth selection in matlab for detail
% One can also manually choose bandwidth, see the following example for
% detail

% default kernel: epanechnikov kernel; other choices: 'normal', 'box', 
% 'triangle'

% generate 100 normal distribution with mean 0 and variance 1
X = normrnd(0, 1, 100, 1);
% plot histogram
histogram(X);
% use the default method (hd method; epanechnikov kernel; default bandwidth)
% to estimate parameters
hd_theta = Continuous_Disparity(X, 'normal');
display(hd_theta)
% Now we specify the input arguments
ned_theta = Continuous_Disparity(X, 'normal', 'disparities', 'ned', ...
            'kernels', 'normal');
display(ned_theta)
% We can even specify the bandwidth (e.g., bandwidth = 1)
ned_theta2 = Continuous_Disparity(X, 'normal', 1, 'disparities', 'ned', ...
            'kernels', 'normal');
display(ned_theta2)



