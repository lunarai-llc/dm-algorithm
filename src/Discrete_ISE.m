% written on 08/16/2018

% this is to compute practical ISE

% the reference is Discrete associated kernels method and extensions



function [ISE] = Discrete_ISE(Y)
    
    a = tabulate(Y);
    YY = a(:, 1);
    
    % Poisson kernel, lambda0 is h, i.e., bandwidth
    % lambda0 = .1;
    % lambda0 = -log((# Y_i = 0)/(sum_i=1^n e^(-Y_i)))
    lambda0 = -log(a(1,2)/(sum(exp(-Y))));
    % cv method:
    % lambda0 = CV_Discrete_kernel(Y);

    g_n = zeros(length(YY), 1);
    for i = 1:length(YY)
        g_n(i) = mean(poisspdf(YY(i), Y+lambda0)); 
    end
    % g_n = g_n/sum(g_n);
    
    ISE = sum((g_n -(.4*poisspdf(YY, .5)+.6*poisspdf(YY, 10)) ).^2 );

end
