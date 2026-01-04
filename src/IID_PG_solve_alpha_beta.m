% written on 02/22/2018

% this is to solve nonlinear equation alpha_j and beta_j for PG model


function G = IID_PG_solve_alpha_beta(x, r_j, Y, nn_j)

G(1) = sum( r_j.*psi(Y+x(1))) - nn_j*(psi(x(1)) - log(x(2)/(x(2)+1)));

G(2) = sum(r_j.*Y*x(2)) - nn_j*x(1);

end


% function F = root2d(x)
% 
% F(1) = exp(-exp(-(x(1)+x(2)))) - x(2)*(1+x(1)^2);
% F(2) = x(1)*cos(x(2)) + x(2)*sin(x(1)) - 0.5;


