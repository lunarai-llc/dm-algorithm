% written on 09/10/2022

% Update on 09/10/2022: This version simplifies the previous version: 
% there is no choice of discrete kernel, only empirical density is selected

% matlab tutorial for creating packages: 
% https://www.mathworks.com/help/matlab/matlab_prog/create-and-share-custom-matlab-toolboxes.html

% This is the discrete disparity method, specifically it implements the
% HD and NED method

% This includes the Binomial, Poisson, Geometric, Negative Binomial, 
% Poisson-Gamma, Poisson-Lognormal model

% The initial estimate is set to me moment estimate

% The optimization method is quasi newton

% Usage: Discrete_Disparity(X, 'hd', 'poiss')

% The standard way to write is shown in Continuous_Disparity, which shows
% how to parse optional arguments.

function [ttheta] = Discrete_Disparity(Y, disparity, model)

argin = inputParser;
argin.addRequired('Y', @isnumeric);
argin.addRequired('disparity', @(x) strcmpi(x,'hd') || ...
    strcmpi(x,'ned'));
argin.addRequired('model', @(x) strcmpi(x,'binom') || ...
    strcmpi(x,'poiss') || strcmpi(x,'geo') || strcmpi(x,'nb') || ...
    strcmpi(x,'pg') || strcmpi(x,'pl') || strcmpi(x,'binom') );

argin.addParameter('aa', 1, @(x) isnumeric(x) && x>0);

% nonparametric part
a = tabulate(Y);
YY = a(:, 1);
g_n = a(:, 3)/100;

options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');  

% disparity method
if strcmp(disparity,'hd')
    if strcmp(model, 'binom')
        fun1 = @(xx) -sum((binopdf(YY, length(YY), xx).*g_n).^.5); 
        x0 = mean(Y)/length(YY);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'poiss')
        fun1 = @(xx) -sum((poisspdf(YY, xx).*g_n).^.5); 
        x0 = mean(Y);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'geo')
        fun1 = @(xx) -sum((geopdf(YY, xx).*g_n).^.5); 
        x0 = 1/(mean(Y)+1);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'nb')
        fun1 = @(xx) -sum((nbinpdf(YY, xx(1), xx(2)).*g_n).^.5); 
        x0 = [mean(Y)/var(Y) mean(Y)^2/var(Y)/(1-mean(Y)/var(Y))];
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'pg')
        fun1 = @(xx) -sum((nbinpdf(YY, xx(1), xx(2)/(xx(2)+1)).*g_n).^.5);
        x0 = [mean(Y)^2/(var(Y)-mean(Y)) mean(Y)/(var(Y)-mean(Y))];
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'pl')
        fun0 = @(lambda, mu, sigma) poisspdf(YY, lambda).* ...
                                    lognpdf(lambda, mu, sigma);
        fun1 = @(xx) -sum((integral(@(lambda)fun0(lambda, xx(1), xx(2)), 0, Inf, 'ArrayValued', true) ...
                      .*g_n).^.5);
        x0 = [log(mean(Y))-log((var(Y)-mean(Y))/mean(Y)^2 +1)/2 log((var(Y)-mean(Y))/mean(Y)^2 +1)];
        ttheta = fminunc(fun1, x0, options);
    else
        error('model not recogonized.');
    end
elseif strcmp(disparity,'ned')
    if strcmp(model, 'binom')
        fun1 = @(xx) sum(exp(-binopdf(YY, length(YY), xx)./g_n).*g_n); 
        x0 = mean(Y)/length(YY);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'poiss')
        fun1 = @(xx) sum(exp(-poisspdf(YY, xx)./g_n).*g_n); 
        x0 = mean(Y);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'geo')
        fun1 = @(xx) sum(exp(-geopdf(YY, xx)./g_n).*g_n); 
        x0 = 1/(mean(Y)+1);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'nb')
        fun1 = @(xx) sum(exp(-nbinpdf(YY, xx(1), xx(2))./g_n).*g_n); 
        x0 = [mean(Y)/var(Y) mean(Y)^2/var(Y)/(1-mean(Y)/var(Y))];
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'pg')
        fun1 = @(xx) sum(exp(-nbinpdf(YY, xx(1), xx(2)/(xx(2)+1))./g_n).*g_n); 
        x0 = [mean(Y)^2/(var(Y)-mean(Y)) mean(Y)/(var(Y)-mean(Y)) ];
        ttheta = fminunc(fun1, x0, options);  
    elseif strcmp(model, 'pl')
        fun0 = @(lambda, mu, sigma) poisspdf(YY, lambda).* ...
                                    lognpdf(lambda, mu, sigma);
        fun1 = @(xx) sum(exp(-integral(@(lambda)fun0(lambda, xx(1), xx(2)), ...
                         0, Inf, 'ArrayValued', true)./g_n).*g_n);
        x0 = [log(mean(Y))-log((var(Y)-mean(Y))/mean(Y)^2 +1)/2 log((var(Y)-mean(Y))/mean(Y)^2 +1)];
        ttheta = fminunc(fun1, x0, options);
    else
        error('model not recogonized.');
    end
else
    error('disparity not recogonized.');
end

end
