% written on 09/10/2022

% parameter structure: pi_1, pi_2, ..., pi_K; theta_1, theta_2, ..., theta_K

% Y is the discrete data

% it can automatically specify the number of parameters

% COMMENT: WORKS WELL

function [z] = Mixture_Model_No_Regress(Y, theta0, disparities, models)

argin = inputParser;

argin.addRequired('Y', @isnumeric);
argin.addRequired('disparities', @(x) strcmpi(x,'em') || strcmpi(x,'hd') || ...
    strcmpi(x,'ned'));
argin.addRequired('models', @(x)  ...
    strcmpi(x,'mixPois') || strcmpi(x,'mixPG') || strcmpi(x,'mixPL'));

n = length(Y);

% specify marginal disparity
a = tabulate(Y);  
% a(any(a==0,2),:) = [];
g_n = a(:, 3)/100; 
YY = a(:, 1);

options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');  

% convergence flag
found = 0;
% convergence tolerance
TOL = 1e-6;
theta = theta0;

iter = 0;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % %           Poisson Mixture   % % % % % % % % % % % % % %      
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

if strcmp(models, 'mixPois')
    % determine K and initialize parameters
    K  = 1/2*(length(theta0));
    pi = theta(1:K);
    lambda = theta((K+1):(2*K));
    dlogl = 1;
    if strcmp(disparities,'em')
    %%%%%%%%%%%%%%%%%%%%%%%%    EM method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        while(found == 0 && iter <= 1000)         
%             pi_old = pi;
%             lambda_old = lambda;        
            dlogl_old = dlogl;
            % update posteria probability:
            bb1 = zeros(n, K);        
            for k = 1:K            
                bb1(:, k) = pi(k)*poisspdf(Y, lambda(k));
            end
            tau = zeros(n, K); denom_bb1 = nansum(bb1, 2);
            for k = 1:K
               tau(:, k) = bb1(:, k)./denom_bb1;
            end
            pi = nansum(tau)/n;
            % update lambdas
            for k = 1:K
                lambda(k) = nansum(Y.*tau(:, k))/nansum(tau(:, k));
            end
            
            %  another method for getting lambda          
            %             for k = 1:K                            
            %                 Logl = @(ttheta) -nanmean(tau(:,k).*log(poisspdf(Y, ttheta)));            
            %                 x1 = lambda(k);    
            %                 lambda(k) = fminunc(Logl, x1, options);                                                    
            %             end                              
            % stopping criteria:
%             error1 = max(abs(pi - pi_old));
%             error2 = max(abs(lambda-lambda_old));        
%             diff1 = max([error1 error2]);   

            f_theta = nansum(bb1);
            
            dlogl = nanmean(log(pi*f_theta'));        
            diff1 = abs(dlogl - dlogl_old);
            
            iter = iter + 1;
            display(iter) 
            display([pi lambda])
            display(diff1)
            if ( diff1 <= TOL)
                found = 1;
            end                      
        end
        % calculate information matrix:
        bb1 = zeros(n, K);        
        for k = 1:K            
            bb1(:, k) = pi(k)*poisspdf(Y, lambda(k));
        end
        tau = zeros(n, K); denom_bb1 = nansum(bb1, 2);
        for k = 1:K
           tau(:, k) = bb1(:, k)./denom_bb1;
        end     
        % variance for pi:
        var_pi = zeros(1, K); 
        for k = 1:K
            var_pi(k) = nansum((tau(:,k)/pi(k)).^2);
        end
        var_pi = 1./var_pi;   
        % variance for lambda    
        var_lambda = zeros(1, K);
        for k = 1:K        
            var_lambda(1, k) = nansum((tau(:,k).*(-1 + Y/lambda(k))).^2);
        end   
        var_lambda = 1./var_lambda;
    %%%%%%%%%%%%%%%%%%%%%%%%    HD method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    elseif strcmp(disparities,'hd')
        dhd = 2;
        while(found == 0 && iter <= 1000)         
%             pi_old = pi;
%             lambda_old = lambda;               
            dhd_old = dhd;
            % update posteria probability:
            bb1 = zeros(length(YY), K);        
            for k = 1:K                
                bb1(:, k) = pi(k)*poisspdf(YY, lambda(k)); 
            end
            tau = zeros(length(YY), K); denom_bb1 = sum(bb1, 2);
            for k = 1:K
               tau(:, k) = bb1(:, k)./denom_bb1;
            end 
            H = zeros(1, K);
            % update lambdas
            for k = 1:K 
                x1 = lambda(k); 
                HD_func = @(ttheta) -nansum((poisspdf(YY, ttheta).*g_n.*tau(:, k)).^.5);                                               
                [lambda(k), H(k)] = fminunc(HD_func, x1, options);                                             
            end  
            pi = H.^2/nansum(H.^2);                        

            % stopping criteria:
%             error1 = max(abs(pi - pi_old));
%             error2 = max(abs(lambda-lambda_old));        
%             diff1 = max([error1 error2]);          

            dhd = -sum(H.*(pi.^.5));       
            diff1 = abs(dhd - dhd_old);
            
            iter = iter + 1;
            display(iter)      
            display([pi lambda])     
            if ( diff1 <= TOL)
                found = 1;
            end                   
        end
        % calculate information matrix:
        bb1 = zeros(length(YY), K);        
        for k = 1:K            
            bb1(:, k) = pi(k)*poisspdf(YY, lambda(k));
        end
        tau = zeros(length(YY), K); denom_bb1 = nansum(bb1, 2);
        for k = 1:K
           tau(:, k) = bb1(:, k)./denom_bb1;
        end         
        W = zeros(length(YY), K); 
        for k = 1:K
            W(:, k) = (pi(k)*poisspdf(YY, lambda(k))./(tau(:,k).*g_n)).^.5;        
        end
        % variance for pi:
        var_pi = zeros(1, K); 
        for k = 1:K
            var_pi(k) = n*nansum((W(:, k).*tau(:,k)/pi(k)).^2.*g_n);
        end
        var_pi = 1./var_pi;
        % variance for lambda
        var_lambda = zeros(1, K);
        for k = 1:K        
            var_lambda(1, k) = n*nansum((W(:, k).*tau(:,k).*(-1+YY/lambda(k))).^2.*g_n);                    
        end   
        var_lambda = 1./var_lambda;
        
    %%%%%%%%%%%%%%%%%%%%%%%%    VNED method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    elseif strcmp(disparities,'vned')
          dvned = 1;
          while(found == 0 && iter <= 1000)         
%               pi_old = pi;
%               lambda_old = lambda;      
              dvned_old = dvned;
              % update posteria probability:
              bb1 = zeros(length(YY), K);        
              for k = 1:K                  
                  bb1(:, k) = pi(k)*poisspdf(YY, lambda(k));                   
              end
              tau = zeros(length(YY), K); denom_bb1 = sum(bb1, 2);
              for k = 1:K
                 tau(:, k) = bb1(:, k)./denom_bb1;
              end 
              % update weight
              W = zeros(length(YY), K); 
              weight_tau = zeros(length(YY), K); 
              for k = 1:K
                  W(:, k) = exp(-pi(k).*poisspdf(YY, lambda(k))./(tau(:,k).*g_n)) ...
                              .*(pi(k).*poisspdf(YY, lambda(k))./(tau(:,k).*g_n)); 
                  weight_tau(:, k) = tau(:, k).*W(:, k);
              end
              % update pi              
              denom_pi = nansum(nansum(weight_tau, 2).*g_n); 
              for k = 1:K
                  pi(k) = nansum(weight_tau(:, k).*g_n)/denom_pi;
              end                               
              % update lambdas
              vned = zeros(1, K);
              for k = 1:K   
                  x1 = lambda(k); 
                  VNED_func = @(ttheta) nansum(exp(-pi(k).*poisspdf(YY, ttheta)./(tau(:, k).*g_n)) ...
                                                   .*tau(:, k).*g_n);                       
                  [lambda(k), vned(k)] = fminunc(VNED_func, x1, options);                   
              end                                      
              
              % stopping criteria:
%               error1 = max(abs(pi - pi_old));
%               error2 = max(abs(lambda-lambda_old));         
%               diff1 = max([error1 error2]);       
                            
              dvned = nansum(vned);        
              diff1 = abs(dvned - dvned_old);
                
              iter = iter + 1;
              display(iter);      
              display([pi lambda]);
              display(diff1);
              
              if (diff1 <= TOL)
                  found = 1;
              end                   
          end
          % calculate information matrix:
            bb1 = zeros(length(YY), K);        
            for k = 1:K            
                bb1(:, k) = pi(k)*poisspdf(YY, lambda(k));
            end
            tau = zeros(length(YY), K); denom_bb1 = nansum(bb1, 2);
            for k = 1:K
               tau(:, k) = bb1(:, k)./denom_bb1;
            end 
            W = zeros(length(YY), K); 
            for k = 1:K
                W(:, k) = exp(-pi(k)*poisspdf(YY, lambda(k))./(tau(:,k).*g_n)).* ...
                             (pi(k)*poisspdf(YY, lambda(k))./(tau(:,k).*g_n));        
            end
            % variance for pi:
            var_pi = zeros(1, K); 
            for k = 1:K
                var_pi(k) = n*nansum((W(:, k).*tau(:,k)/pi(k)).^2.*g_n);
            end
            var_pi = 1./var_pi;
            % variance for lambdas    
            var_lambda = zeros(1, K);
            for k = 1:K        
                var_lambda(1, k) = n*nansum((W(:, k).*tau(:,k).*(-1+YY/lambda(k))).^2.*g_n);   
            end   
            var_lambda = 1./var_lambda;
    else
        error('disparity not recogonized.');        
    end
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % %     Poisson Gamma Mixture   % % % % % % % % % % % % % %      
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    

elseif strcmp(models, 'mixPG')
    % determine K and initialize parameters
    K  = 1/3*(length(theta0));
    pi = theta(1:K);
    alpha = theta((K+1):(2*K));
    beta = theta((2*K+1):(3*K));
    
    if strcmp(disparities,'em')
    %%%%%%%%%%%%%%%%%%%%%%%%    EM method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        dlogl = 1;
        while(found == 0 && iter <= 1000)         
%             pi_old = pi;
%             alpha_old = alpha;
%             beta_old = beta;
            
            dlogl_old = dlogl;
                 
            % update posteria probability:
            bb1 = zeros(n, K);        
            for k = 1:K            
                bb1(:, k) = pi(k)*nbinpdf(Y, alpha(k), beta(k)/(beta(k)+1));
            end
            tau = zeros(n, K); denom_bb1 = nansum(bb1, 2);
            for k = 1:K
               tau(:, k) = bb1(:, k)./denom_bb1;
            end
            sum_tau = nansum(tau);  % sum_tau is a vector
            pi = sum_tau/n;
            
            % solve alpha_k and beta_k together.
            for k = 1:K
                fun = @(x) IID_PG_solve_alpha_beta(x, tau(:, k), Y, sum_tau(k));        
                para = fsolve(fun, [alpha(k) beta(k)]);        
                alpha(k) = para(1);
                beta(k)  = para(2);        
            end   
            
            % stopping criteria:
%             error1 = max(abs(pi - pi_old));
%             error2 = max(abs(alpha-alpha_old));        
%             error3 = max(abs(beta-beta_old));   
%             
%             diff1 = max([error1 error2 error3]);    

            f_theta = nansum(bb1);
            
            dlogl = nanmean(log(pi*f_theta'));        
            diff1 = abs(dlogl - dlogl_old);
            
            iter = iter + 1;
            display(iter) 
            display([pi alpha beta])     
            if ( diff1 <= TOL)
                found = 1;
            end                      
        end
        % calculate information matrix:
        bb1 = zeros(n, K);        
        for k = 1:K            
            bb1(:, k) = pi(k)*nbinpdf(Y, alpha(k), beta(k)/(beta(k)+1));
        end
        tau = zeros(n, K); denom_bb1 = nansum(bb1, 2);
        for k = 1:K
           tau(:, k) = bb1(:, k)./denom_bb1;
        end     
        
        % variance for pi:
        var_pi = zeros(1, K); 
        for k = 1:K
            var_pi(k) = nansum((tau(:,k)/pi(k)).^2);
        end
        var_pi = 1./var_pi;           
        % variance for alphas and betas    
        var_alpha = zeros(1, K);
        var_beta = zeros(1, K);
        for k = 1:K    
            var_alpha(1, k) = nansum((tau(:, k).*(psi(Y+alpha(k)) - psi(alpha(k)) + log(beta(k)/(beta(k)+1)))).^2);
            var_beta(1, k) = nansum( (tau(:, k).*(Y*beta(k) - alpha(k)/(beta(k)*(beta(k)+1))) ).^2 );
        end   
        var_alpha = 1./var_alpha;
        var_beta = 1./var_beta;
    %%%%%%%%%%%%%%%%%%%%%%%%    HD method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    elseif strcmp(disparities,'hd')
        dhd = 2;
        while(found == 0 && iter <= 1000)         
%             pi_old = pi;
%             alpha_old = alpha;
%             beta_old = beta;      
            
            dhd_old = dhd;
            % update posteria probability:
            bb1 = zeros(length(YY), K);        
            for k = 1:K                
                bb1(:, k) = pi(k)*nbinpdf(YY, alpha(k), beta(k)/(beta(k)+1));                
            end
            tau = zeros(length(YY), K); denom_bb1 = sum(bb1, 2);
            for k = 1:K
               tau(:, k) = bb1(:, k)./denom_bb1;
            end                         
            H = zeros(1, K);
            % update alphas and betas
            for k = 1:K 
                H1 = @(ttheta) -nansum((nbinpdf(YY, ttheta(1), ttheta(2)/(ttheta(2)+1)) ...
                                        .*g_n.*tau(:,k)).^.5 ); 
                x1 = [alpha(k) beta(k)];    
                [phi1, H(k)] = fminunc(H1, x1, options);    
                alpha(k) = phi1(1);
                beta(k)  = phi1(2);                                                                     
            end  
            
            pi = H.^2/nansum(H.^2);

            % stopping criteria:
%             error1 = max(abs(pi - pi_old));
%             error2 = max(abs(alpha-alpha_old));        
%             error3 = max(abs(beta-beta_old));   
%             diff1 = max([error1 error2 error3]); 
            
            dhd = -sum(H.*(pi.^.5));  
            diff1 = max(abs(dhd - dhd_old));                         
            
            iter = iter + 1;
            display(iter);      
            display([pi alpha beta]);
            display(diff1);
            
            if ( diff1 <= TOL)
                found = 1;
            end                   
        end
        % calculate information matrix:
        bb1 = zeros(length(YY), K);        
        for k = 1:K            
            bb1(:, k) = pi(k)*nbinpdf(YY, alpha(k), beta(k)/(beta(k)+1));                
        end
        tau = zeros(length(YY), K); denom_bb1 = nansum(bb1, 2);
        for k = 1:K
           tau(:, k) = bb1(:, k)./denom_bb1;
        end             
        
        W = zeros(length(YY), K); 
        for k = 1:K
            W(:, k) = (pi(k)*nbinpdf(YY, alpha(k), beta(k)/(beta(k)+1))./(tau(:,k).*g_n)).^.5;        
        end             
        
        % variance for pi:
        var_pi = zeros(1, K); 
        for k = 1:K
            var_pi(k) = n*nansum((W(:, k).*tau(:,k)/pi(k)).^2.*g_n);
        end
        var_pi = 1./var_pi;
        
        % variance for alpha and beta
        var_alpha = zeros(1, K);
        var_beta = zeros(1, K);
        for k = 1:K        
            var_alpha(1, k) = n* nansum((W(:, k).*tau(:, k).*(psi(YY+alpha(k)) ...
                                        - psi(alpha(k)) + log(beta(k)/(beta(k)+1)))).^2.*g_n);
            var_beta(1, k) = n* nansum((W(:, k).*tau(:, k).* ...
                                       (YY*beta(k)-alpha(k)/(beta(k)*(beta(k)+1))) ).^2.*g_n);
        end   
        var_alpha = 1./var_alpha;
        var_beta = 1./var_beta;
        
    %%%%%%%%%%%%%%%%%%%%%%%%    VNED method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    elseif strcmp(disparities,'vned')
          dvned = 1;
          while(found == 0 && iter <= 1000)         
%               pi_old = pi;
%               alpha_old = alpha;
%               beta_old = beta;         
              dvned_old = dvned;
              % update posteria probability:
              bb1 = zeros(length(YY), K);        
              for k = 1:K
                  bb1(:, k) = pi(k)*nbinpdf(YY, alpha(k), beta(k)/(beta(k)+1));                
              end
              tau = zeros(length(YY), K); denom_bb1 = sum(bb1, 2);
              for k = 1:K
                 tau(:, k) = bb1(:, k)./denom_bb1;
              end 
              
              % update weight
              W = zeros(length(YY), K); 
              weight_tau = zeros(length(YY), K); 
              for k = 1:K
                  W(:, k) = exp(-pi(k)*nbinpdf(YY, alpha(k), beta(k)/(beta(k)+1))./(tau(:,k).*g_n)).* ...
                             (pi(k)*nbinpdf(YY, alpha(k), beta(k)/(beta(k)+1))./(tau(:,k).*g_n)); 
                  weight_tau(:, k) = tau(:, k).*W(:, k);
              end
              % update pi              
              denom_pi = nansum(nansum(weight_tau, 2).*g_n); 
              for k = 1:K
                  pi(k) = nansum(weight_tau(:, k).*g_n)/denom_pi;
              end                               
              % update alphas and betas
              vned = zeros(1, K);
              for k = 1:K   
                  x1 = [alpha(k) beta(k)];    
                  VNED_func = @(ttheta) nansum(exp(-pi(k)*nbinpdf(YY, ttheta(1), ttheta(2)/(ttheta(2)+1))./(tau(:, k).*g_n)) ...
                                                .*tau(:, k).*g_n);                                               
                  [phi1, vned(k)] = fminunc(VNED_func, x1, options);    
                  alpha(k) = phi1(1);
                  beta(k)  = phi1(2);                                                                                       
              end                                    
              
              % stopping criteria:
%               error1 = max(abs(pi - pi_old));
%               error2 = max(abs(alpha-alpha_old));        
%               error3 = max(abs(beta-beta_old));   
%             
%               diff1 = max([error1 error2 error3]);         

              dvned = nansum(vned);
              diff1 = abs(dvned - dvned_old);
              
              iter = iter + 1;
              display(iter);      
              display([pi alpha beta]);
              display(diff1);
              if (diff1 <= TOL)
                  found = 1;
              end                   
          end
          % calculate information matrix:
            bb1 = zeros(length(YY), K);        
            for k = 1:K            
                bb1(:, k) = pi(k)*nbinpdf(YY, alpha(k), beta(k)/(beta(k)+1));  
            end
            tau = zeros(length(YY), K); denom_bb1 = nansum(bb1, 2);
            for k = 1:K
               tau(:, k) = bb1(:, k)./denom_bb1;
            end 
            W = zeros(length(YY), K); 
            for k = 1:K
                W(:, k) = exp(-pi(k)*nbinpdf(YY, alpha(k), beta(k)/(beta(k)+1))./(tau(:,k).*g_n)).* ...
                             (pi(k)*nbinpdf(YY, alpha(k), beta(k)/(beta(k)+1))./(tau(:,k).*g_n));        
            end
            % variance for pi:
            var_pi = zeros(1, K); 
            for k = 1:K
                var_pi(k) = n*nansum((W(:, k).*tau(:,k)/pi(k)).^2.*g_n);
            end
            var_pi = 1./var_pi;
             
            % variance for alphas and betas
            var_alpha = zeros(1, K);
            var_beta = zeros(1, K);
            for k = 1:K        
                var_alpha(1, k) = n* nansum((W(:, k).*tau(:, k).*(psi(YY+alpha(k)) ...
                                            - psi(alpha(k)) + log(beta(k)/(beta(k)+1)))).^2.*g_n);
                var_beta(1, k) = n* nansum((W(:, k).*tau(:, k).* ...
                                           (YY*beta(k)-alpha(k)/(beta(k)*(beta(k)+1))) ).^2.*g_n);
            end   
            var_alpha = 1./var_alpha;
            var_beta = 1./var_beta;
    else
        error('disparity not recogonized.');        
    end    
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % %     Poisson LogNormal Mixture   % % % % % % % % % % % % % %      
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    

elseif strcmp(models, 'mixPL')
    % determine K and initialize parameters
    K  = 1/3*(length(theta0));
    pi = theta(1:K);
    mu = theta((K+1):(2*K));
    sigma = theta((2*K+1):(3*K));        
    
    if strcmp(disparities,'em')
    %%%%%%%%%%%%%%%%%%%%%%%%    EM method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        % Density function without integral
        ffun1 = @(uu, phi) poisspdf(Y, uu).*lognpdf(uu, phi(1), phi(2)); 
        % First and derivatives of the density w.r.t. mu and sigma
        ffun_mu = @(uu, phi) poisspdf(Y, uu).*lognpdf(uu, phi(1), phi(2)) ...
                                .*(log(uu)-phi(1))/(phi(2))^2; 
        ffun_sigma = @(uu, phi) poisspdf(Y, uu).*lognpdf(uu, phi(1), phi(2)) ...
                                .*(-1/phi(2) + (log(uu) -phi(1))^2/(phi(2))^3 ); 
                            
        dlogl = 1;
        while(found == 0 && iter <= 1000)         
%             pi_old = pi;
%             mu_old = mu;
%             sigma_old = sigma;
            
            dlogl_old = dlogl;
                 
            % update posteria probability:
            bb1 = zeros(n, K);        
            for k = 1:K            
                bb1(:, k) = pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                           0, 100, 'ArrayValued', true);     
            end
            tau = zeros(n, K); denom_bb1 = nansum(bb1, 2);
            for k = 1:K
               tau(:, k) = bb1(:, k)./denom_bb1;
            end
            sum_tau = nansum(tau);  % sum_tau is a vector
            pi = sum_tau/n;
            
            % solve mu_k and sigma_k together.
            f_theta = zeros(n, k);
            for k = 1:K                
                Logl = @(xx) -nanmean(tau(:, k).*log(integral(@(uu)ffun1(uu, xx), ...
                                                  0, 100, 'ArrayValued', true) ));  
                x1 = [mu(k) sigma(k)];    
                [phi1, ~] = fminunc(Logl, x1, options);    
                mu(k) = phi1(1); 
                sigma(k) = phi1(2);                
                f_theta(:, k) = integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                           0, 100, 'ArrayValued', true);                                                           
            end   
            
            % stopping criteria:
%             error1 = max(abs(pi - pi_old));
%             error2 = max(abs(mu-mu_old));        
%             error3 = max(abs(sigma-sigma_old));               
%             diff1 = max([error1 error2 error3]);  
            
            % stopping criteria:
            dlogl = nanmean(log(pi*f_theta'));        
            diff1 = abs(dlogl - dlogl_old);
            
            iter = iter + 1;
            display(iter); 
            display([pi mu sigma]);
            display(diff1);
            if ( diff1 <= TOL)
                found = 1;
            end                      
        end
        % calculate information matrix:
        bb1 = zeros(n, K);        
        for k = 1:K            
            bb1(:, k) = pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                           0, 100, 'ArrayValued', true); 
        end
        tau = zeros(n, K); denom_bb1 = nansum(bb1, 2);
        for k = 1:K
           tau(:, k) = bb1(:, k)./denom_bb1;
        end     
        
        % variance for pi:
        var_pi = zeros(1, K); 
        for k = 1:K
            var_pi(k) = nansum((tau(:,k)/pi(k)).^2);
        end
        var_pi = 1./var_pi;           
        % variance for mus and sigmas    
        var_mu = zeros(1, K);
        var_sigma = zeros(1, K);
        for k = 1:K    
            var_mu(1, k) = nansum((tau(:, k).*(integral(@(uu)ffun_mu(uu, [mu(k) sigma(k)]), ...
                                               0, 100, 'ArrayValued', true) ...
                                               ./integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                               0, 100, 'ArrayValued', true) )).^2);
            var_sigma(1, k) = nansum((tau(:, k).*(integral(@(uu)ffun_sigma(uu, [mu(k) sigma(k)]), ...
                                                  0, 100, 'ArrayValued', true) ...
                                                  ./integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                                  0, 100, 'ArrayValued', true) )).^2);
        end   
        var_mu = 1./var_mu;
        var_sigma = 1./var_sigma;
    %%%%%%%%%%%%%%%%%%%%%%%%    HD method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    elseif strcmp(disparities,'hd')
        % Density function without integral
        ffun1 = @(uu, phi) poisspdf(YY, uu).*lognpdf(uu, phi(1), phi(2)); 
        % First and derivatives of the density w.r.t. mu and sigma
        ffun_mu = @(uu, phi) poisspdf(YY, uu).*lognpdf(uu, phi(1), phi(2)) ...
                                .*(log(uu)-phi(1))/(phi(2))^2; 
        ffun_sigma = @(uu, phi) poisspdf(YY, uu).*lognpdf(uu, phi(1), phi(2)) ...
                                .*(-1/phi(2) + (log(uu) -phi(1))^2/(phi(2))^3 ); 
        dhd = 2;
        
        while(found == 0 && iter <= 1000)         
%             pi_old = pi;
%             mu_old = mu;
%             sigma_old = sigma;  
            
            dhd_old = dhd;
            % update posteria probability:
            bb1 = zeros(length(YY), K);        
            for k = 1:K                
                bb1(:, k) = pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                           0, 100, 'ArrayValued', true);               
            end
            tau = zeros(length(YY), K); denom_bb1 = sum(bb1, 2);
            for k = 1:K
               tau(:, k) = bb1(:, k)./denom_bb1;
            end                         
            H = zeros(1, K);
            % update mus and sigmas
            for k = 1:K 
                H1 = @(ttheta) -nansum((integral(@(uu)ffun1(uu, ttheta), ...
                                        0, 100, 'ArrayValued', true) ...
                                        .*g_n.*tau(:,k)).^.5 ); 
                x1 = [mu(k) sigma(k)];    
                [phi1, H(k)] = fminunc(H1, x1, options);    
                mu(k) = phi1(1);
                sigma(k)  = phi1(2);                                                                     
            end  
            
            pi = H.^2/nansum(H.^2);

            % stopping criteria:
%             error1 = max(abs(pi - pi_old));
%             error2 = max(abs(mu-mu_old));        
%             error3 = max(abs(sigma-sigma_old));   
%             
%             diff1 = max([error1 error2 error3]); 
            
            dhd = -sum(H.*(pi.^.5));         
            % stopping criteria:
            diff1 = abs(dhd - dhd_old);

            iter = iter + 1;
            display(iter);      
            display([pi mu sigma]);
            display(diff1);
            
            if ( diff1 <= TOL)
                found = 1;
            end                   
        end
        % calculate information matrix:
        bb1 = zeros(length(YY), K);        
        for k = 1:K            
            bb1(:, k) = pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                           0, 100, 'ArrayValued', true);                 
        end
        tau = zeros(length(YY), K); denom_bb1 = nansum(bb1, 2);
        for k = 1:K
           tau(:, k) = bb1(:, k)./denom_bb1;
        end             
        
        W = zeros(length(YY), K); 
        for k = 1:K
            W(:, k) = pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                           0, 100, 'ArrayValued', true);         
        end             
        
        % variance for pi:
        var_pi = zeros(1, K); 
        for k = 1:K
            var_pi(k) = n*nansum((W(:, k).*tau(:,k)/pi(k)).^2.*g_n);
        end
        var_pi = 1./var_pi;
        
        % variance for mu and sigma
        var_mu = zeros(1, K);
        var_sigma = zeros(1, K);
        for k = 1:K        
            var_mu(1, k) = n* nansum((W(:, k).*tau(:, k).*(integral(@(uu)ffun_mu(uu, [mu(k) sigma(k)]), ...
                                               0, 100, 'ArrayValued', true) ...
                                               ./integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                               0, 100, 'ArrayValued', true))).^2.*g_n);
            var_sigma(1, k) = n* nansum((W(:, k).*tau(:, k).* ...
                                               (integral(@(uu)ffun_sigma(uu, [mu(k) sigma(k)]), ...
                                               0, 100, 'ArrayValued', true) ...
                                               ./integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                               0, 100, 'ArrayValued', true)) ).^2.*g_n);
        end   
        var_mu = 1./var_mu;
        var_sigma = 1./var_sigma;
        
    %%%%%%%%%%%%%%%%%%%%%%%%    VNED method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    elseif strcmp(disparities,'vned')
          % Density function without integral
          ffun1 = @(uu, phi) poisspdf(YY, uu).*lognpdf(uu, phi(1), phi(2)); 
          % First and derivatives of the density w.r.t. mu and sigma
          ffun_mu = @(uu, phi) poisspdf(YY, uu).*lognpdf(uu, phi(1), phi(2)) ...
                               .*(log(uu)-phi(1))/(phi(2))^2; 
          ffun_sigma = @(uu, phi) poisspdf(YY, uu).*lognpdf(uu, phi(1), phi(2)) ...
                                  .*(-1/phi(2) + (log(uu) -phi(1))^2/(phi(2))^3 );   
                              
          dvned = 1;
          while(found == 0 && iter <= 1000)  
%               pi_old = pi;
%               mu_old = mu;
%               sigma_old = sigma;  
              dvned_old = dvned;
              % update posteria probability:
              bb1 = zeros(length(YY), K);        
              for k = 1:K
                  bb1(:, k) = pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                             0, 100, 'ArrayValued', true);                
              end
              tau = zeros(length(YY), K); denom_bb1 = sum(bb1, 2);
              for k = 1:K
                 tau(:, k) = bb1(:, k)./denom_bb1;
              end 
              
              % update weight
              W = zeros(length(YY), K); 
              weight_tau = zeros(length(YY), K); 
              for k = 1:K
                  W(:, k) = exp(-pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                       0, 100, 'ArrayValued', true)./(tau(:,k).*g_n)).* ...
                               (pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                0, 100, 'ArrayValued', true)./(tau(:,k).*g_n)); 
                  weight_tau(:, k) = tau(:, k).*W(:, k);
              end
              % update pi              
              denom_pi = nansum(nansum(weight_tau, 2).*g_n); 
              for k = 1:K
                  pi(k) = nansum(weight_tau(:, k).*g_n)/denom_pi;
              end                               
              % update alphas and betas
              
              vned = zeros(1, K);
              for k = 1:K   
                  x1 = [mu(k) sigma(k)];    
                  VNED_func = @(ttheta) nansum(exp(-pi(k)*integral(@(uu)ffun1(uu, ttheta), ...
                                               0, 100, 'ArrayValued', true)./(tau(:, k).*g_n)) ...
                                               .*tau(:, k).*g_n);                                               
                  [phi1, vned(k)] = fminunc(VNED_func, x1, options);    
                  mu(k) = phi1(1);
                  sigma(k)  = phi1(2);                                     
              end                                                                
              
              % stopping criteria:
%               error1 = max(abs(pi - pi_old));
%               error2 = max(abs(mu-mu_old));        
%               error3 = max(abs(sigma-sigma_old));   
% 
%               diff1 = max([error1 error2 error3]); 
                
              dvned = nansum(vned);

              diff1 = max(abs(dvned - dvned_old));  

              iter = iter + 1;
              display(iter);      
              display([pi mu sigma]);
              display(diff1);
            
              if (diff1 <= TOL)
                  found = 1;
              end                   
          end
          % calculate information matrix:
            bb1 = zeros(length(YY), K);        
            for k = 1:K            
                bb1(:, k) = pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                           0, 100, 'ArrayValued', true); 
            end
            tau = zeros(length(YY), K); denom_bb1 = nansum(bb1, 2);
            for k = 1:K
               tau(:, k) = bb1(:, k)./denom_bb1;
            end 
            W = zeros(length(YY), K); 
            for k = 1:K
                W(:, k) =exp(-pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                    0, 100, 'ArrayValued', true)./(tau(:,k).*g_n)) ...
                         .*(pi(k)*integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                            0, 100, 'ArrayValued', true)./(tau(:,k).*g_n));        
            end
            % variance for pi:
            var_pi = zeros(1, K); 
            for k = 1:K
                var_pi(k) = n*nansum((W(:, k).*tau(:,k)/pi(k)).^2.*g_n);
            end
            var_pi = 1./var_pi;            
            % variance for mu and sigma
            var_mu = zeros(1, K);
            var_sigma = zeros(1, K);
            for k = 1:K        
                var_mu(1, k) = n* nansum((W(:, k).*tau(:, k).*(integral(@(uu)ffun_mu(uu, [mu(k) sigma(k)]), ...
                                               0, 100, 'ArrayValued', true) ...
                                               ./integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                               0, 100, 'ArrayValued', true))).^2.*g_n);
                var_sigma(1, k) = n* nansum((W(:, k).*tau(:, k).* ...
                                               (integral(@(uu)ffun_sigma(uu, [mu(k) sigma(k)]), ...
                                               0, 100, 'ArrayValued', true) ...
                                               ./integral(@(uu)ffun1(uu, [mu(k) sigma(k)]), ...
                                               0, 100, 'ArrayValued', true)) ).^2.*g_n);
            end   
            var_mu = 1./var_mu;
            var_sigma = 1./var_sigma;
    else
        error('disparity not recogonized.');        
    end       
    
end


% Output:

if strcmp(models, 'mixPois')
    z.pi = pi;
    z.theta = [pi lambda(:)'];

    z.var_pi = var_pi;
    z.var_lambda = var_lambda;
    
elseif strcmp(models, 'mixPG')
    z.pi = pi;
    z.theta = [pi alpha beta];

    z.var_pi = var_pi;
    z.var_alpha = var_alpha;
    z.var_beta = var_beta;
    
elseif strcmp(models, 'mixPL')
    z.pi = pi;
    z.theta = [pi mu sigma];

    z.var_pi = var_pi;
    z.var_mu = var_mu;
    z.var_sigma = var_sigma;    
    
end
    

end







