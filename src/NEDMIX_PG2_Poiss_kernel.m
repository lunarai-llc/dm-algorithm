% written on 08/20/2018

% true_pi = 0.3; true_alpha1 = 10; true_beta1 = 1;
% true_alpha2 = 1; true_beta2 = 2;
% 
% n = 5000; 
% 
% true_theta = [true_pi true_alpha1 true_beta1 true_alpha2 true_beta2];
% theta0 = true_theta;
% Y = Generate_Y_iid_PG(true_theta, n);
% % initial value:
% theta0 = INI_PG(Y);
% 
% % initial value:
% theta0 = INI_PG(F);


function [z] = NEDMIX_PG2_Poiss_kernel(Y, theta0) 

z = [];

% k = 2;

% lambda0 = .2;

a = tabulate(Y);
YY = a(:, 1);
% cc = a(:, 2);

% g_n = a(:, 3)/100;

g_n = Poisson_kernel(Y);

iter = 0;

pi0 = theta0(1); alpha1 = theta0(2); beta1 = theta0(3);
                 alpha2 = theta0(4); beta2 = theta0(5);      
 

rr1 = pi0*nbinpdf(YY, alpha1, beta1/(beta1+1)) ...      
      +(1-pi0)*nbinpdf(YY, alpha2, beta2/(beta2+1));

r(:, 1) = pi0*nbinpdf(YY, alpha1, beta1/(beta1+1))./rr1;

r(:, 2) = (1-pi0)*nbinpdf(YY, alpha2, beta2/(beta2+1))./rr1; 
                                      
ned = sum(exp(-pi0*nbinpdf(YY, alpha1, beta1/(beta1+1))./(g_n.*r(:, 1))) ...
               .*(g_n.*r(:, 1))) + ...
      sum(exp(-(1-pi0)*nbinpdf(YY, alpha2, beta2/(beta2+1))./(g_n.*r(:, 2))) ...
                .*(g_n.*r(:, 2)));

theta = theta0;   

% convergence flag
found = 0;

% convergence tolerance
TOL = 1e-6;

pi    = theta(1);
alpha = [theta(2) theta(4)]';
beta  = [theta(3) theta(5)]';

while( found == 0 )
    
    ned_old = ned;
    
    % save current beta
    % theta_old = theta;
    % logl_old  = logl;        
       
    %==============%
    %    E Step    % 
    %==============%
    
    r = zeros(length(YY), 2);        % since k = 2 it is two component mixture
    
    r(:, 1) = pi*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1))./(pi*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1)) ...
              +(1-pi)*nbinpdf(YY, alpha(2), beta(2)/(beta(2)+1)) );

    r(:, 2) = (1-pi)*nbinpdf(YY, alpha(2), beta(2)/(beta(2)+1))./(pi*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1)) ...
              +(1-pi)*nbinpdf(YY, alpha(2), beta(2)/(beta(2)+1)) ); 
          
    % pi = sum(cc.*r(:, 1)/n);      
    
    %==============%
    %    M Step    % 
    %==============%        
               
    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');           
    
    NED1 = @(x) sum(exp(-(pi*nbinpdf(YY, x(1), x(2)/(x(2)+1))./(g_n.*r(:, 1)))).*g_n.*r(:, 1) );     
    
    % this is another type of objective function.
    %     NED1 = @(x) sum(exp(-((g_n.*r(:, 1))./(pi*nbinpdf(YY, x(1), x(2)/(x(2)+1))) )) ...
    %                         .*(pi*nbinpdf(YY, x(1), x(2)/(x(2)+1)))  );  
    
    x1 = [alpha(1) beta(1)];    
    [phi1, ~] = fminunc(NED1, x1, options);    
    alpha(1) = phi1(1);
    beta(1)  = phi1(2);
    % hd1 is the hellinger distance.
    % hd1 = -sum((nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1)).*g_n.*r(:,1)).^.5 );
    
    %============================%
    %========  PI1  works =======%
    %============================%
    
    %     PI1 = @(aa) sum(exp(-((aa.*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1))+ ...
    %                     (1-aa).*nbinpdf(YY, alpha(2), beta(2)/(beta(2)+1)) )./g_n)) ...
    %                     .*g_n );                     
    % 
    %     pi =  fminunc(PI1, pi, options);  
            
            
    NED2 = @(xx) sum(exp(-((1-pi)*nbinpdf(YY, xx(1), xx(2)/(xx(2)+1))./(g_n.*r(:, 2)))).*g_n.*r(:, 2) );                                                      
    
    x2 = [alpha(2) beta(2)];    
    [phi2, ~] = fminunc(NED2, x2, options);    
    alpha(2) = phi2(1);
    beta(2)  = phi2(2); 
    
    
    pi = sum(exp(-pi*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1))./(g_n.*r(:, 1))) ...
                 .*(pi*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1))))/  ...
         (sum(exp(-pi*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1))./(g_n.*r(:, 1))) ...
                 .*(pi*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1)))) + ...
          sum(exp(-(1-pi)*nbinpdf(YY, alpha(2), beta(2)/(beta(2)+1))./(g_n.*r(:, 2))) ...
                 .*((1-pi)*nbinpdf(YY, alpha(2), beta(2)/(beta(2)+1))))  );
             
    
    ned = sum(exp(-pi*nbinpdf(YY, alpha(1), beta(1)/(beta(1)+1))./(g_n.*r(:, 1))) ...
               .*(g_n.*r(:, 1))) + ...
          sum(exp(-(1-pi)*nbinpdf(YY, alpha(2), beta(2)/(beta(2)+1))./(g_n.*r(:, 2))) ...
                .*(g_n.*r(:, 2)));
        
    theta = [pi alpha(1) beta(1) alpha(2) beta(2)]; 
    
    iter = iter + 1;
    
    display(iter)
    
    if ( max(abs(ned - ned_old)) <= TOL )
        found = 1;
    end
    
    %     if ( max(abs(logl - logl_old)) <= TOL )
    %         found = 1;
    %     end
    
    display(theta)
    
end

z.theta = theta;

z.mixture = 2;

z.ned = ned;

end

