% written on 08/20/2018

% this is to use EM type algorithm for HD method to implement iid PG model

% true_pi = 0.3; true_alpha1 = 10; true_beta1 = 1;
% true_alpha2 = 1; true_beta2 = 2;
% 
% n = 500; 
% 
% true_theta = [true_pi true_alpha1 true_beta1 true_alpha2 true_beta2];
% theta0 = true_theta;
% Y = Generate_Y_iid_PG(true_theta, n);

% EM algorithm 

% COMMENT: WORKS WELL, it is essentially the same as HMIX algorithm

% COMMENT: THIS IS Example 1 kernel


function [z] = HMIX_PG2_Ex1_kernel(Y, theta0, aa) 

z = [];

% k = 2;

% lambda0 = .5;

a = tabulate(Y);
YY = a(:, 1);
% cc = a(:, 2);

% % empirical density
% g_n = a(:, 3)/100;
% aa = 1;
g_n = Example1_kernel(Y, aa);

iter = 0;

pi0 = theta0(1); alpha1 = theta0(2); beta1 = theta0(3);
alpha2 = theta0(4); beta2 = theta0(5);

rr1 = (pi0*nbinpdf(YY, alpha1, beta1/(beta1+1)) ...
          +(1-pi0)*nbinpdf(YY, alpha2, beta2/(beta2+1)) );

r(:, 1) = pi0*nbinpdf(YY, alpha1, beta1/(beta1+1))./rr1;

r(:, 2) = (1-pi0)*nbinpdf(YY, alpha2, beta2/(beta2+1))./rr1;

h = -(pi0^.5*sum((nbinpdf(YY, alpha1, beta1/(beta1+1)).*r(:, 1).*g_n).^.5) ...
      +(1-pi0)^.5*sum((nbinpdf(YY, alpha2, beta2/(beta2+1)).*r(:, 2).*g_n).^.5));
       
theta = theta0;   

% convergence flag
found = 0;

% convergence tolerance
TOL = 10e-7;

pi    = theta(1);
alpha = [theta(2) theta(4)]';
beta  = [theta(3) theta(5)]'; 

while( found == 0 )
        
    h_old = h;    
    
    %==============%
    %    E Step    % 
    %==============%
    
    r = zeros(length(YY), 2);        % since k = 2 it is two component mixture
    
    rr1 = (pi*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1)) ...
           +(1-pi)*nbinpdf(YY, alpha(2), beta(2)/(beta(2)+1)) );
       
    r(:, 1) = pi*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1))./rr1;

    r(:, 2) = (1-pi)*nbinpdf(YY, alpha(2), beta(2)/(beta(2)+1))./rr1;                                            
                              
    %==============%
    %    M Step    % 
    %==============%        
               
    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');           
    
    H1 = @(x) -sum((nbinpdf(YY, x(1), x(2)/(x(2)+1)).*g_n.*r(:, 1)).^.5); 
    x1 = [alpha(1) beta(1)];    
    [phi1, f1] = fminunc(H1, x1, options);    
    alpha(1) = phi1(1);
    beta(1)  = phi1(2);

    H2 = @(xx) -sum((nbinpdf(YY, xx(1), xx(2)/(xx(2)+1)).*g_n.*r(:, 2)).^.5);      
    x2 = [alpha(2) beta(2)];    
    [phi2, f2] = fminunc(H2, x2, options);    
    alpha(2) = phi2(1);
    beta(2)  = phi2(2); 
    
    pi = f1^2/(f1^2+f2^2); 
    
    h = -(pi^.5*sum((nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1)).*r(:, 1).*g_n).^.5) ...
          +(1-pi)^.5*sum((nbinpdf(YY, alpha(2), beta(2)/(beta(2)+1)).*r(:, 2).*g_n).^.5));            
    
    theta = [pi alpha(1) beta(1) alpha(2) beta(2)]; 
    
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

