% written on 08/16/2018

% this is the poisson kernel

% Discrete associated kernels method and extensions

% COMMENT: this code is correct

function [g_n] = Poisson_kernel(Y)

a = tabulate(Y);
YY = a(:, 1);

% Poisson kernel, lambda0 is h, i.e., bandwidth
% lambda0 = .1;
% lambda0 = -log((# Y_i = 0)/(sum_i=1^n e^(-Y_i)))

lambda0 = -log(a(1,2)/(sum(exp(-Y))));

% cv method:
% lambda0 = CV_Poisson_kernel(Y);

% lambda0 = 0.05;

g_n = zeros(length(YY), 1);
for i = 1:length(YY)
    g_n(i) = mean(poisspdf(Y, YY(i)+lambda0)); 
end

% nn = 20;
% g_n = zeros(nn, 1);
% for i = 1:nn
%     g_n(i) = mean(poisspdf(Y, i-1+lambda0)); 
% end
% 
% g_n = g_n(YY+1)/sum(g_n);

end


