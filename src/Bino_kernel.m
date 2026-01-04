% written on 08/17/2018

% this is the Binomial discrete kernel

% reference: Discrete associated kernels method and extensions

% COMMENT: this code is correct

function [g_n] = Bino_kernel(Y)

a = tabulate(Y);
YY = a(:, 1);
% cc = a(:, 2);

%==============================================%
%   this is one criterion to choose bandwidth  %
%==============================================%

ffun = @(x) sum(((1-x)./(Y+1)).^(Y+1)) - a(1, 2) ;
options = optimoptions('fsolve','Display','off');
x0 = .5;
lambda0 = fsolve(ffun,x0,options);

% lambda0 = 0.05;

%==============================================%
%   % % cv method to choose bandwidth          %
%==============================================%
% lambda0 = CV_Bino_kernel(Y);
% lambda0 = .1;

% g_n = zeros(length(YY), 1);
% for i = 1:length(YY)    
%     g_n(i) = mean(binopdf(YY(i), Y+1, (Y+lambda0)./(Y+1)));                
% end

% g_n = zeros(length(YY), 1);
% for i = 1:length(YY)
%     cc = 0;
%     for j = 1:length(Y)
%         if YY(i) - (Y(j)+1) <=0
%             cc = cc + binopdf(YY(i), Y(j)+1, (Y(j)+lambda0)/(Y(j)+1));            
%         end
%     end
%     g_n(i) = cc/length(Y);
% end

% g_n = g_n/sum(g_n);

g_n = zeros(length(YY), 1);
for i = 1:length(YY)    
    g_n(i) = mean(binopdf(Y, YY(i)+1, (YY(i)+lambda0)./(YY(i)+1)));                
end

% nn = 20;
% g_n = zeros(nn, 1);
% for i = 1:nn   
%     g_n(i) = mean(binopdf(Y, i-1+1, (i-1+lambda0)./(i-1+1)));                
% end
% 
% g_n = g_n(YY+1)/sum(g_n);


end




