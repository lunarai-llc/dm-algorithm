%% ISE with different kernels for .4*Poi(.5) + .6*Poi(10)
% change sample size n to get different results and then make the plots
rng(666);

MM = 5000;

n = 100;

theta0 = [.4 .5 10];

lambda0 = 0.1; % lambda0 is the bandwidth

ise_empirical = zeros(MM, 1);

ise_example1_aa1  = zeros(MM, 1);  
ise_example1_aa2  = zeros(MM, 1);  

ise_bino = zeros(MM, 1);

ise_poisson = zeros(MM, 1);

ise_nb     = zeros(MM, 1);

tic
parfor kk = 1:MM
    
    Y = Generate_Poi2(theta0, n);
                
    ise_example1_aa1(kk) = ISE_two_Poisson_Example1(Y, 1, lambda0);
    ise_example1_aa2(kk) = ISE_two_Poisson_Example1(Y, 2, lambda0);
    ise_bino(kk) = ISE_two_Poisson_Bino_kernel(Y, lambda0);
    ise_poisson(kk) = ISE_two_Poisson_Poisson_kernel(Y, lambda0);
    ise_nb(kk) = ISE_two_Poisson_NB_kernel(Y, lambda0);
    ise_empirical(kk) = ISE_two_Poisson_Empirical(Y);

end
toc

display([mean(ise_empirical)*1000 std(ise_empirical)*1000])
display([mean(ise_example1_aa1)*1000 std(ise_example1_aa1)*1000])
display([mean(ise_example1_aa2)*1000 std(ise_example1_aa2)*1000])
display([mean(ise_bino)*1000 std(ise_bino)*1000])
display([mean(ise_poisson)*1000 std(ise_poisson)*1000])
display([mean(ise_nb)*1000 std(ise_nb)*1000])

