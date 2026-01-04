% written on 08/17/2018

% this is example 1 ISE

% Specifically, ISE = sum_{x in N} (g_n(x) f(x; theta) )^2

% Reference: Discrete associated kernels method and extensions

% COMMENT: this code is correct

function [ISE] = ISE_two_Poisson_Example1(Y, aa, lambda0)
    
    %     a = tabulate(Y);
    %     YY = a(:, 1);
    %     g_nn = a(:, 3)/100;

    %============================================%
    %   take all possible values into account    %
    %============================================%

    %     aa = 1;
    nn = 30;
    
    % % cv method:
    % lambda0 = CV_Example1_kernel(Y, aa);
    
    
    % lambda0 = .1;

    % aa is a in the example
    % pp is P(a, h) in the example
    % lambda0 is h in the example
    % bb = sum_k=0^a k^lambda0
    % lambda0 = 0.1;
    
    if nargin < 2 || isempty(lambda0)

        % Option 2: 
        lambda0 = tri_zero_mass_bandwidth(Y, aa, false);
    end    

    bb = 0;
    for k = 0:aa
        bb = bb + k^lambda0;
    end
    pp = (2*aa + 1)*(aa + 1)^lambda0 - 2*bb;

    % g_n = zeros(length(YY), 1);
    % for i = 1:length(YY)
    %     g_n(i) = mean(((aa+1)^lambda0 - ((Y-YY(i))<=aa)* ((Y-YY(i))>=0 )*(abs(Y-YY(i))).^lambda0 )./pp );
    % end

    g_n = zeros(nn, 1);
    for i = 1:nn
        cc = 0;
        for j = 1:length(Y)
            if  abs(Y(j) - (i-1)) <= aa
                cc = cc + ((aa+1)^lambda0 - (abs(Y(j)-(i-1)))^lambda0 )/pp;
            end
        end
        g_n(i) = cc/length(Y);
    end
    
    xx = 0:(nn-1);
    
%      fun = poisspdf(xx, 2);
    fun = .4*poisspdf(xx, .5) + .6*poisspdf(xx, 10);
% fun = .3*nbinpdf(xx, 10, 1/(1+1)) + .7*nbinpdf(xx, 1, 2/(2+1));

    ISE = sum((g_n -  fun').^2);

    %============================================%
    %  only take the sample value into account   %
    %============================================%

%     a = tabulate(Y);
%     YY = a(:, 1);
%     
%     g_n = Example1_kernel(Y);
%     fun = poisspdf(YY, 2);
%     
%     ISE = sum((g_n -  fun).^2);

end