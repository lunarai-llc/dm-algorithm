% written on 08/17/2018

% this is ISE for Poisson kernel.

% Specifically, ISE = sum_{x in N} (g_n(x) f(x; theta) )^2

% Reference: Discrete associated kernels method and extensions

function [ISE] = ISE_two_Poisson_Poisson_kernel(Y, lambda0)


a = tabulate(Y);

% If lambda0 is not provided, set a default
    if nargin < 2 || isempty(lambda0)
        % Option 1: use CV_Bino_kernel
        % lambda0 = CV_Bino_kernel(Y);

        % Option 2: 
       lambda0 = -log(a(1,2)/(sum(exp(-Y))));
    end


nn = 30;

% lambda0 = 0.1;

% lambda0 = CV_Poisson_kernel(Y);
% lambda0 = 0.7;

xx = 0:(nn-1);

g_n = zeros(nn, 1);
for x = 1:nn
    g_n(x) = mean(poisspdf(Y, x-1+lambda0)); 
end

fun = .4*poisspdf(xx, .5) + .6*poisspdf(xx, 10);

% fun = poisspdf(xx, 2);

% fun = .3*nbinpdf(xx, 10, 1/(1+1)) + .7*nbinpdf(xx, 1, 2/(2+1));

% a = tabulate(Y);
% YY = a(:, 1);
% 
% lambda0 = 0.1;
% 
% % xx = 0:(nn-1);
% 
% g_n = zeros(length(YY), 1);
% for i = 1:nn
%     g_n(i) = mean(poisspdf(YY(i), Y+lambda0)); 
% end

% g_n = Poisson_kernel(Y);

% fun = true_pdf(Y);

ISE = sum((g_n -  fun').^2);
    
end