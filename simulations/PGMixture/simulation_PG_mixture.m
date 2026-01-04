% written on 08/20/2018



%% two component PG mixture

rng(666)

theta0 = [.3 10 1 1 2];

MM = 5000;

n = 1000;

true_theta = repmat(theta0, MM, 1);

em_theta = zeros(MM, 5);

hd_ex1   = zeros(MM, 5);
hd_ex2   = zeros(MM, 5);
hd_bino  = zeros(MM, 5);
hd_poiss = zeros(MM, 5);
hd_nb    = zeros(MM, 5);
hd_emp   = zeros(MM, 5);

ned_ex1   = zeros(MM, 5);
ned_ex2   = zeros(MM, 5);
ned_bino  = zeros(MM, 5);
ned_poiss = zeros(MM, 5);
ned_nb    = zeros(MM, 5);
ned_emp   = zeros(MM, 5);

tic;

parfor k = 1:MM
    
      Y = Generate_Y_iid_PG2(theta0, n);
      
      ini_theta = INI_PG2(Y);
    
      em_theta(k, :) = EM_iid_PG2(Y, ini_theta);
    
      % example 1 kernel
      c1 = HMIX_PG2_Ex1_kernel(Y, ini_theta, 1);
      hd_ex1(k, :) = c1.theta;
    
      c2 = HMIX_PG2_Ex1_kernel(Y, ini_theta, 2);
      hd_ex2(k, :) = c2.theta;    
    
      c3 = HMIX_PG2_Bino_kernel(Y, ini_theta);
      hd_bino(k, :) = c3.theta;
    
      c4 = HMIX_PG2_Poiss_kernel(Y, ini_theta);
      hd_poiss(k, :) = c4.theta;
    
      c5 = HMIX_PG2_NB_kernel(Y, ini_theta);
      hd_nb(k, :) = c5.theta;
    
      c6 = EM_iid_PG2_HD(Y, ini_theta);
      hd_emp(k, :) = c6.theta;
    
    
      d1 =  NEDMIX_PG2_Ex1_kernel(Y, ini_theta, 1);   
      ned_ex1(k, :) = d1.theta;
      
      d2 =  NEDMIX_PG2_Ex1_kernel(Y, ini_theta, 2);   
      ned_ex2(k, :) = d2.theta;

      d3 = NEDMIX_PG2_Bino_kernel(Y, ini_theta);  
      ned_bino(k, :) = d3.theta;

      d4 = NEDMIX_PG2_Poiss_kernel(Y, ini_theta);   
      ned_poiss(k, :) = d4.theta;
    
      d5 = NEDMIX_PG2_NB_kernel(Y, ini_theta);  
      ned_nb(k, :) = d5.theta;

      d6 = EM_iid_PG2_NED(Y, ini_theta);
      ned_emp(k, :) = d6.theta;
     
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

