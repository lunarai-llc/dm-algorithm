% written on 08/16/2018

% this is the kernel proposed by Kokonendji

% This is the example 1 from the reference below

% Reference: Discrete associated kernels method and extensions

% COMMENT: this code is correct

function [g_n] = Example1_kernel(Y, aa)

a = tabulate(Y);
YY = a(:, 1);
% cc = a(:, 2);

% % cv method:
% lambda0 = CV_Example1_kernel(Y, aa);
% lambda0 = .05;

lambda0 = tri_zero_mass_bandwidth(Y, aa, false);

% aa is a in the example
% pp is P(a, h) in the example
% lambda0 is h in the example
% bb = sum_k=0^a k^lambda0
% lambda0 = 0.1;

bb = 0;
for k = 0:aa
    bb = bb + k^lambda0;
end
pp = (2*aa + 1)*(aa + 1)^lambda0 - 2*bb;

% g_n = zeros(length(YY), 1);
% for i = 1:length(YY)
%     g_n(i) = mean(((aa+1)^lambda0 - ((Y-YY(i))<=aa)* ((Y-YY(i))>=0 )*(abs(Y-YY(i))).^lambda0 )./pp );
% end

g_n = zeros(length(YY), 1);
for i = 1:length(YY)
    cc = 0;
    for j = 1:length(Y)
        if  abs(Y(j) - YY(i)) <= aa
            cc = cc + ((aa+1)^lambda0 - (abs(Y(j)-YY(i)))^lambda0 )/pp;
        end
    end
    g_n(i) = cc/length(Y);
end

% g_n = g_n/sum(g_n);

end




