% written on 08/16/2018

% this is the negative binomial kernel

% Reference: Discrete associated kernels method and extensions

function [g_n] = NB_kernel(Y)

a = tabulate(Y);
YY = a(:, 1);
% cc = a(:, 2);

% negative binomial kernel
ffun = @(x) sum(((Y+1)./(2*Y+1+x)).^(Y+1)) - a(1, 2) ;
options = optimoptions('fsolve','Display','off');
x0 = .5;
lambda0 = fsolve(ffun,x0,options);

% lambda0 = 0.05;

% % cv method:
% lambda0 = CV_Discrete_kernel(Y);

% lambda0 = .95;
% [x,fval,exitflag,output] = fsolve(ffun,x0,options);

g_n = zeros(length(YY), 1);
for i = 1:length(YY)
    g_n(i) = mean(nbinpdf(Y, YY(i)+1, (YY(i)+1)./(2*YY(i)+1+lambda0)));
end

% nn = 20;
% g_n = zeros(nn, 1);
% for i = 1:nn
%     g_n(i) = mean(nbinpdf(Y, i-1+1, (i-1+1)./(2*(i-1)+1+lambda0)));
% end
% g_n = g_n(YY+1)/sum(g_n);

end
