% written on 08/14/2018

% this is to use HMIX for two component Poisson mixture

% this is to use discrete poisson kernel

% another reference: Discrete associated kernels method and extensions

% COMMENT: WORKS WELL, not very sensitive to initial values

function [z] = HMIX_two_Poisson_Poisson_kernel(Y, theta0)

z = [];

a = tabulate(Y);
YY = a(:, 1);
% cc = a(:, 2);

% binomial kernel
% g_n = Bino_kernel(Y);

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

iter = 0;

% pi0 = theta0(1); lambda1 = theta0(2); lambda2 = theta0(3);
% 
% r1 = pi0*poisspdf(YY, lambda1);
% r2 = (1-pi0)*poisspdf(YY, lambda2);
% 
% r(:, 1) = r1./(r1+r2);
% r(:, 2) = r2./(r1+r2);
% 
% h = -(pi0^.5*sum((poisspdf(YY, lambda1).*r(:, 1).*g_n).^.5) ...
%       +(1-pi0)^.5*sum((poisspdf(YY, lambda2).*r(:, 2).*g_n).^.5));

h = 2;
       
theta = theta0;   

% convergence flag
found = 0;

% convergence tolerance
TOL = 10e-7;

pi    = theta(1);
lambda = [theta(2) theta(3)]';

while( found == 0 )
        
    h_old = h;    
    
    %==============%
    %    E Step    % 
    %==============%
    
    r = zeros(length(YY), 2);        % since k = 2 it is two component mixture
    
    r1 = pi*poisspdf(YY, lambda(1));
    r2 = (1-pi)*poisspdf(YY, lambda(2));
    
    r(:, 1) = r1./(r1+r2);
    r(:, 2) = r2./(r1+r2);                                                   
                              
    %==============%
    %    M Step    % 
    %==============%        
               
    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');           
    
    H1 = @(x) -sum((poisspdf(YY, x(1)).*g_n.*r(:, 1)).^.5); 
    x1 = lambda(1);   
    [phi1, f1] = fminunc(H1, x1, options);    
    lambda(1) = phi1(1);    

    H2 = @(xx) -sum((poisspdf(YY, xx(1)).*g_n.*r(:, 2)).^.5);      
    x2 = lambda(2);    
    [phi2, f2] = fminunc(H2, x2, options);    
    lambda(2) = phi2(1);
    
    % update pi
    pi = f1^2/(f1^2+f2^2); 
    
    h = -(pi^.5*sum((poisspdf(YY, lambda(1)).*r(:, 1).*g_n).^.5) ...
          +(1-pi)^.5*sum((poisspdf(YY, lambda(2)).*r(:, 2).*g_n).^.5));            
    
    theta = [pi lambda(1) lambda(2)]; 
    
    iter = iter + 1;
    
    display(iter)        
    
    if ( max(abs(h - h_old)) <= TOL )
        found = 1;
    end               
    
    display(theta)
    
end

z.theta = theta;

z.mixture = 2;

z.hd = 2+2*h;

end

