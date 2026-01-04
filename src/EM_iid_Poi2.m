% written on 04/20/2018

% this is EM poisson algorithm

function [ttheta] = EM_iid_Poi2(Y, theta0) 

% k = 2;

iter = 0;

n = length(Y);

lambda = zeros(1, 2);

pi0 = theta0(1); lambda(1) = theta0(2); lambda(2) = theta0(3);

logl = mean(log(pi0*poisspdf(Y,lambda(1)) ...
           +(1-pi0)*poisspdf(Y,lambda(2)) ));       
       
theta = theta0;   

% convergence flag
found = 0;

% convergence tolerance
TOL = 1e-5;

r = zeros(n, 2);

while( found == 0 )
    
    % save current beta
    % theta_old = theta;
    logl_old  = logl;
    
    pi    = theta(1);
    lambda = [theta(2) theta(3)];
       
    %==============%
    %    E Step    % 
    %==============%
    
    % r = zeros(n, 2);        % since k = 2 it is two component mixture
    
    r(:, 1) = pi*poisspdf(Y,lambda(1))./(pi*poisspdf(Y,lambda(1)) ...
              +(1-pi)*poisspdf(Y,lambda(2)) );

    r(:, 2) = (1-pi)*poisspdf(Y,lambda(2))./(pi*poisspdf(Y,lambda(1)) ...
              +(1-pi)*poisspdf(Y,lambda(2)) );

    nn = sum(r);
    
    %==============%
    %    M Step    % 
    %==============%
    
    % update proportion pi
    pi = nn/n;
    
    % solve lambda_j.
    lambda(1) = sum(Y.*r(:, 1))/sum(r(:, 1));
    lambda(2) = sum(Y.*r(:, 2))/sum(r(:, 2));
                 
    logl = mean(log(pi(1)*poisspdf(Y,lambda(1)) ...
                +(1-pi(1))*poisspdf(Y,lambda(2)) ));
    
    theta = [pi(1) lambda(1) lambda(2)]; 
    
    iter = iter + 1;
    
    display(iter)
    
    if ( max(abs(logl - logl_old)) <= TOL )
        found = 1;
    end
    
    display(theta)
    
end

ttheta = theta;

end
