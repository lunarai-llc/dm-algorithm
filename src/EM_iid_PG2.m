% written on 02/22/2018


% this is to use EM algorithm for IID PG model

% this is using solving the nonlinear equation to solve the problem

% generate IID PG random variables

% true_pi = 0.3; true_alpha1 = 10; true_beta1 = 1;
% true_alpha2 = 1; true_beta2 = 2;
% 
% n = 10000; 
% 
% true_theta = [true_pi true_alpha1 true_beta1 true_alpha2 true_beta2];
% theta0 = true_theta;
% Y = Generate_Y_iid_PG(true_theta, n);

% 
% % initial value:
% theta0 = INI_PG(F);

% EM algorithm 

% COMMENT: solve alpha and beta together works well. also solve beta first
% then update alpha also works well.

% WORKS WELL

function [ttheta] = EM_iid_PG2(Y, theta0) 

k = 2;

iter = 0;

n = length(Y);

pi0 = theta0(1); alpha1 = theta0(2); beta1 = theta0(3);
alpha2 = theta0(4); beta2 = theta0(5);

logl = mean(log(pi0*nbinpdf(Y, alpha1, beta1/(beta1+1)) ...
           +(1-pi0)*nbinpdf(Y, alpha2, beta2/(beta2+1)) ));       
       
theta = theta0;   

% convergence flag
found = 0;

% convergence tolerance
TOL = 1e-6;

while( found == 0 )
    
    % save current beta
    % theta_old = theta;
    logl_old  = logl;
    
    pi    = theta(1);
    alpha = [theta(2) theta(4)]';
    beta  = [theta(3) theta(5)]';
       
    %==============%
    %    E Step    % 
    %==============%
    
    r = zeros(n, 2);        % since k = 2 it is two component mixture
    
    r(:, 1) = pi*nbinpdf(Y, alpha(1), beta(1)/(beta(1)+1))./(pi*nbinpdf(Y, alpha(1), beta(1)/(beta(1)+1)) ...
              +(1-pi)*nbinpdf(Y, alpha(2), beta(2)/(beta(2)+1))  );

    r(:, 2) = (1-pi)*nbinpdf(Y, alpha(2), beta(2)/(beta(2)+1))./(pi*nbinpdf(Y, alpha(1), beta(1)/(beta(1)+1)) ...
              +(1-pi)*nbinpdf(Y, alpha(2), beta(2)/(beta(2)+1)) );

    nn = sum(r);
    
    %==============%
    %    M Step    % 
    %==============%
    
    % update proportion pi
    pi = nn/n;
    
    % solve alpha_j and beta_j together.
    for j = 1:k
        fun = @(x) IID_PG_solve_alpha_beta(x, r(:, j), Y, nn(j));        
        para = fsolve(fun, [alpha(j) beta(j)]);        
        alpha(j) = para(1);
        beta(j)  = para(2);        
    end
    
    %     % solve for alpha_j first then solve beta_j
    %     % COMMENT: WORKS WELL
    %     for j = 1:k
    %         fun = @(x) IID_PG_solve_alpha(x, r(:, j), Y, nn(j));
    %         alpha(j) = fsolve(fun, alpha(j)); 
    %     end
    % 
    %     for j = 1:k       
    %         beta(j) = nn(j)*alpha(j)/sum(r(:, j).*Y); 
    %     end
    
    
    % solve beta_j first then solve alpha_j
    % this method is questionable, but the simulation works well.
    
    %     % solve for beta_j
    %     for j = 1:k
    %         beta(j) = nn(j)*alpha(j)/sum(r(:, j).*Y);    
    %     end
    % 
    %     % solve for alpha_j   
    %     % the reference is: https://www.mathworks.com/help/optim/ug/fsolve.html#butbmfz-5
    %     for j = 1:k
    %         fun = @(x) IID_PG_solve_alpha(x, r(:, j), Y, nn(j), beta(j));
    %         alpha(j) = fsolve(fun, alpha(j));  
    %     end
    
    logl = mean(log(pi(1)*nbinpdf(Y, alpha(1), beta(1)/(beta(1)+1)) ...
                +(1-pi(1))*nbinpdf(Y, alpha(2), beta(2)/(beta(2)+1)) ));
    
    theta = [pi(1) alpha(1) beta(1) alpha(2) beta(2)]; 
    
    iter = iter + 1;
    
    display(iter)
    
    if ( max(abs(logl - logl_old)) <= TOL )
        found = 1;
    end
    
    display(theta)
    
end

ttheta = theta;

end

