% written on 08/17/2018

% this is ISE for empirical value


function [ISE] = ISE_two_Poisson_Empirical(Y)
    
    %============================================%
    %   take all possible values into account    %
    %============================================%

    %     a = tabulate(Y); 
    %     YY = a(:, 1);            
    
    nn = 30;
    g_n = zeros(nn, 1);
    for x = 1:nn
        g_n(x) = mean(Y == (x-1));
    end
      
    xx = 0:(nn-1);
    
%      fun = poisspdf(xx, 2);
    fun = .4*poisspdf(xx, .5) + .6*poisspdf(xx, 10);
%     fun = .3*nbinpdf(xx, 10, 1/(1+1)) + .7*nbinpdf(xx, 1, 2/(2+1));

    ISE = sum((g_n -  fun').^2);

    %============================================%
    %  only take the sample value into account   %
    %============================================%

%     a = tabulate(Y);
%     YY = a(:, 1);
%     g_n = a(:, 3)/100;     
%     fun = poisspdf(YY, 2);
%     
%     ISE = sum((g_n -  fun).^2);

end