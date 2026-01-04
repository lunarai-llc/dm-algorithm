% written on 08/17/2018

% this is ISE for Binomial kernel.

% Specifically, ISE = sum_{x in N} (g_n(x) f(x; theta) )^2

% Reference: Discrete associated kernels method and extensions

function [ISE] = ISE_two_Poisson_Bino_kernel(Y, lambda0)

% discrete kernel
% lambda0 = CV_Bino_kernel(Y);
% lambda0 = 0.7;

a = tabulate(Y);

% If lambda0 is not provided, set a default
    if nargin < 2 || isempty(lambda0)
        % Option 1: use CV_Bino_kernel
        % lambda0 = CV_Bino_kernel(Y);

        % Option 2: 
        ffun = @(x) sum(((1-x)./(Y+1)).^(Y+1)) - a(1, 2) ;
        options = optimoptions('fsolve','Display','off');
        x0 = .5;
        lambda0 = fsolve(ffun,x0,options);
    end

xx = 0:(max(Y)+1);

g_n = zeros(length(xx), 1);
for x = 1:(max(Y)+2)            
    g_n(x) = mean(binopdf(Y, x-1+1, (x-1+lambda0)/(x-1+1)));                     
end

% g_n = zeros(max(Y)+1, 1);
% for x = 1:(max(Y)+1)
%     g_n(x) = mean(binopdf(x-1, Y+1, (Y+lambda0)./(Y+1)));
% end

% fun = poisspdf(xx, 2);

fun = .4*poisspdf(xx, .5) + .6*poisspdf(xx, 10);
% fun = .3*nbinpdf(xx, 10, 1/(1+1)) + .7*nbinpdf(xx, 1, 2/(2+1));

ISE = sum((g_n -  fun').^2);
    
end

