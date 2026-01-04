% written on 08/19/2018

% this is NEDMIX for discrete kernel

% this is two component Poisson

% COMMENT: WORKS WELL

function [ttheta] = NEDMIX_two_Poisson_Poisson_kernel(Y, theta0) 

iter = 0;

% z = [];

a = tabulate(Y);
YY = a(:, 1);
% cc = a(:, 2);

% % empirical density
% g_n = a(:, 3)/100;

% Poisson kernel:
g_n = Poisson_kernel(Y);

% Negative binomial kernel:
% g_n = NB_kernel(Y);

% example 1:
% aa = 1;
% g_n = Example1_kernel(Y, aa);

% example 4:
% g_n = Example4_kernel(Y);

% convergence flag
found = 0;
     
theta = theta0;  

% convergence tolerance
TOL = 10e-7;

pi    = theta(1);
lambda = [theta(2) theta(3)]';

ned = 2; 

while( found == 0 )
        
    ned_old = ned;           
                                                                                    
    %==============%
    %    M Step    % 
    %==============% 
    
    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');    
    
    ned1 = @(x) sum(exp(-pi*poisspdf(YY, x) ...
             ./(pi*poisspdf(YY, lambda(1))./ ...
            (pi*poisspdf(YY, lambda(1))+(1-pi)*poisspdf(YY, lambda(2))).*g_n) ) ...
             .*(pi*poisspdf(YY, lambda(1))./ ...
            (pi*poisspdf(YY, lambda(1))+(1-pi)*poisspdf(YY, lambda(2))).*g_n)) ; 
            
    x1 = lambda(1);
    [phi1, ned1] = fminunc(ned1, x1, options);        
        
    
    ned2 = @(x) sum(exp(-(1-pi)*poisspdf(YY, x) ...
             ./((1-pi)*poisspdf(YY, lambda(2))./ ...
            (pi*poisspdf(YY, lambda(1))+(1-pi)*poisspdf(YY, lambda(2))).*g_n) ) ...
             .*((1-pi)*poisspdf(YY, lambda(2))./ ...
            (pi*poisspdf(YY, lambda(1))+(1-pi)*poisspdf(YY, lambda(2))).*g_n)) ; 
                                                                        
    x2 = lambda(2);    
    [phi2, ned2] = fminunc(ned2, x2, options);        
        
    lambda(1) = phi1(1);  lambda(2) = phi2(1);         
    
    %=============%
    %  update pi  % 
    %=============%
    
    nedd1 = sum(exp(-pi*poisspdf(YY, lambda(1)) ...
             ./(pi*poisspdf(YY, lambda(1))./ ...
            (pi*poisspdf(YY, lambda(1))+(1-pi)*poisspdf(YY, lambda(2))).*g_n) ) ...
             .* pi.*poisspdf(YY, lambda(1)));
                 
         
    nedd2 =  sum(exp(-(1-pi)*poisspdf(YY, lambda(2)) ...
             ./((1-pi)*poisspdf(YY, lambda(2))./ ...
             (pi*poisspdf(YY, lambda(1))+(1-pi)*poisspdf(YY, lambda(2))).*g_n) ) ...
             .*(1-pi).*poisspdf(YY, lambda(2)));                        
    
    pi  = nedd1/(nedd1 + nedd2);
        
    ned = ned1 + ned2;
            
    theta = [pi lambda(1) lambda(2)]; 
    
    iter = iter + 1;
    
    display(iter)
        
    if ( max(abs(ned - ned_old)) <= TOL )
        found = 1;
    end                 
    
    display(theta)
    
end

ttheta = theta;

end







