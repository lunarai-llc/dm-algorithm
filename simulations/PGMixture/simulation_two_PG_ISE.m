


%% all different kernels for .3*PG(10, 1) + .7*PG(1, 2), bandwidth needs to be specified

MM = 5000;

n = 1000;

theta0 = [.3 10 1 1 2];

lambda0 = 0.1; % lambda0 is the bandwidth

ise_empirical = zeros(MM, 1);

ise_example1_aa1  = zeros(MM, 1);  
ise_example1_aa2  = zeros(MM, 1);  

ise_bino = zeros(MM, 1);

ise_poisson = zeros(MM, 1);

ise_nb     = zeros(MM, 1);

tic
parfor kk = 1:MM
    
    Y = Generate_Y_iid_PG2(theta0, n);
                
    ise_example1_aa1(kk) = ISE_two_PG_Example1(Y, 1, lambda0);
    ise_example1_aa2(kk) = ISE_two_PG_Example1(Y, 2, lambda0);
    ise_bino(kk) = ISE_two_PG_Bino_kernel(Y, lambda0);
    ise_poisson(kk) = ISE_two_PG_Poisson_kernel(Y, lambda0);
    ise_nb(kk) = ISE_two_PG_NB_kernel(Y, lambda0);
    ise_empirical(kk) = ISE_two_PG_Empirical(Y);

end
toc


%% display results
display([mean(ise_empirical)*1000 std(ise_empirical)*1000])
display([mean(ise_example1_aa1)*1000 std(ise_example1_aa1)*1000])
display([mean(ise_example1_aa2)*1000 std(ise_example1_aa2)*1000])
display([mean(ise_bino)*1000 std(ise_bino)*1000])
display([mean(ise_poisson)*1000 std(ise_poisson)*1000])
display([mean(ise_nb)*1000 std(ise_nb)*1000])











