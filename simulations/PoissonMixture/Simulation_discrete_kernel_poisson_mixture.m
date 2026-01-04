% written on 08/14/2018

% this is the simulation for Poisson mixture

% Specifically, we compare the empirical distribution with discrete kernel


%% Poisson mixture

rng(666)

theta0 = [.4 .5 10];

MM = 5000;

n = 200;

true_theta = repmat(theta0, MM, 1);

em_theta = zeros(MM, 3);

hd_ex1   = zeros(MM, 3);
hd_ex2   = zeros(MM, 3);
hd_bino  = zeros(MM, 3);
hd_poiss = zeros(MM, 3);
hd_nb    = zeros(MM, 3);
hd_emp   = zeros(MM, 3);

ned_ex1   = zeros(MM, 3);
ned_ex2   = zeros(MM, 3);
ned_bino  = zeros(MM, 3);
ned_poiss = zeros(MM, 3);
ned_nb    = zeros(MM, 3);
ned_emp   = zeros(MM, 3);

tic;

parfor k = 1:MM
    
    Y = Generate_Poi2(theta0, n);
    
    ini_theta = INI_Poi2(Y);
    
    em_theta(k, :) = EM_iid_Poi2(Y, ini_theta);
    
    c1 = HMIX_two_Poisson_Example1_kernel(Y, ini_theta, 1);
    hd_ex1(k, :) = c1.theta;
    
    c2 = HMIX_two_Poisson_Example1_kernel(Y, ini_theta, 2);
    hd_ex2(k, :) = c2.theta;    
    
    c3 = HMIX_two_Poisson_Bino_kernel(Y, ini_theta);
    hd_bino(k, :) = c3.theta;
    
    c4 = HMIX_two_Poisson_Poisson_kernel(Y, ini_theta);
    hd_poiss(k, :) = c4.theta;
    
    c5 = HMIX_two_Poisson_NB_kernel(Y, ini_theta);
    hd_nb(k, :) = c5.theta;
    
    c6 = HMIX_two_Poisson_Empirical_kernel(Y, ini_theta);
    hd_emp(k, :) = c6.theta;
    
    
    ned_ex1(k, :) = NEDMIX_two_Poisson_Example1_kernel(Y, ini_theta, 1);    
    
    ned_ex2(k, :) = NEDMIX_two_Poisson_Example1_kernel(Y, ini_theta, 2);    
    
    ned_bino(k, :) = NEDMIX_two_Poisson_Bino_kernel(Y, ini_theta);
    
    ned_poiss(k, :) = NEDMIX_two_Poisson_Poisson_kernel(Y, ini_theta);    
    
    ned_nb(k, :) = NEDMIX_two_Poisson_NB_kernel(Y, ini_theta);    
    
    ned_emp(k, :) = NEDMIX_two_Poisson_Empirical__kernel(Y, ini_theta);
     
end

toc


%%
ave_em = mean(em_theta);         std_em = std(em_theta);
ave_hd_emp = mean(hd_emp);       std_hd_emp = std(hd_emp);
ave_hd_ex1 = mean(hd_ex1);       std_hd_ex1 = std(hd_ex1);
ave_hd_ex2 = mean(hd_ex2);       std_hd_ex2 = std(hd_ex2);
ave_hd_bino = mean(hd_bino);     std_hd_bino = std(hd_bino);
ave_hd_poiss = mean(hd_poiss);   std_hd_poiss = std(hd_poiss);
ave_hd_nb = mean(hd_nb);         std_hd_nb = std(hd_nb);

ave_ned_emp = mean(ned_emp);       std_ned_emp = std(ned_emp);
ave_ned_ex1 = mean(ned_ex1);       std_ned_ex1 = std(ned_ex1);
ave_ned_ex2 = mean(ned_ex2);       std_ned_ex2 = std(ned_ex2);
ave_ned_bino = mean(ned_bino);     std_ned_bino = std(ned_bino);
ave_ned_poiss = mean(ned_poiss);   std_ned_poiss = std(ned_poiss);
ave_ned_nb = mean(ned_nb);         std_ned_nb = std(ned_nb);


display([ave_em std_em])

display([ave_hd_emp std_hd_emp])
display([ave_hd_ex1 std_hd_ex1])
display([ave_hd_ex2 std_hd_ex2])
display([ave_hd_bino std_hd_bino])
display([ave_hd_poiss std_hd_poiss])
display([ave_hd_nb std_hd_nb])

display([ave_ned_emp std_ned_emp])
display([ave_ned_ex1 std_ned_ex1])
display([ave_ned_ex2 std_ned_ex2])
display([ave_ned_bino std_ned_bino])
display([ave_ned_poiss std_ned_poiss])
display([ave_ned_nb std_ned_nb])









