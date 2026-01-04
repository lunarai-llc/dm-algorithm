% written on 08/17/2018

% this is ISE for Poisson kernel.

% Specifically, ISE = sum_{x in N} (g_n(x) f(x; theta) )^2

% Reference: Discrete associated kernels method and extensions

function [ISE] = ISE_two_Poisson_NB_kernel(Y, lambda0)

nn = 30;

a = tabulate(Y);

% If lambda0 is not provided, set a default
if nargin < 2 || isempty(lambda0)
    % Option 1: use CV_Bino_kernel
    % lambda0 = CV_Bino_kernel(Y);

    % Option 2: 
    % negative binomial kernel
    ffun = @(x) sum(((Y+1)./(2*Y+1+x)).^(Y+1)) - a(1, 2) ;
    options = optimoptions('fsolve','Display','off');
    x0 = .5;
    lambda0 = fsolve(ffun,x0,options);
end


xx = 0:(nn-1);

g_n = zeros(nn, 1);
for x = 1:nn
    g_n(x) = mean(nbinpdf(Y, x-1+1, (x-1+1)./(2*(x-1)+1+lambda0)));
end

% fun = poisspdf(xx, 2);
fun = .4*poisspdf(xx, .5) + .6*poisspdf(xx, 10);
% fun = .3*nbinpdf(xx, 10, 1/(1+1)) + .7*nbinpdf(xx, 1, 2/(2+1));

ISE = sum((g_n -  fun').^2);
    
end



