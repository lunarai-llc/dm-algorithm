% written on 08/21/2018




%% bandwidth h = 0.1 for PG model


k1 = [20 50 100 200 500 1000 2000];

emp = [36.52 14.20 7.383 3.843 1.468 0.774 0.364];
ex1 = [30.46 13.43 8.006 5.061 3.524 2.828 2.507];
ex2 = [30.46 17.30 12.96 10.50 9.552 8.905 8.603];
bino = [18.52 8.806 5.773 4.073 3.220 2.803 2.606];
poiss = [25.28 17.46 15.17 13.87 13.22 12.90 12.77];
nb = [40.68 32.69 30.48 29.29 28.49 28.24 28.15];

%

plot(k1, emp, '-or', 'LineWidth', 1.2)
hold on
plot(k1, ex1, '-*b', 'LineWidth', 1.2)
plot(k1, ex2, '-*g', 'LineWidth', 1.2)
plot(k1, bino, '-*c', 'LineWidth', 1.2)
plot(k1, poiss, '--k', 'LineWidth', 1.2)
plot(k1, nb, '-.b', 'LineWidth', 1.2)


legend('Empirical', 'Triangle a=1', 'Triangle a = 2', 'Binomial', 'Poisson', 'Negative Binomial')

% title('Mean of ISE for PG model h = 0.1')





% %%
% 
% % legend('pnc', 'ave of CSA')
% % ylim([0 .845])
% ylabel('Estimates')
% xlabel('1: \pi_1, 2: \alpha_1, 3: \beta_1, 4: \alpha_2, 5: \beta_2')
% XTick = k1;
% 
% % xticks([1 2 3 4 5])
% % xticklabels({'\pi1','\alpha1','\beta_1', '\alpha_2','\beta_2'})
% set(gca, 'xtick', XTick)
% 
% % for kk=1:numel(k1)
% %       text(k1(kk),ave_csa(kk),['(' num2str(k1(kk)) ',' num2str(ave_csa(kk)) ')'])
% % end




%% bandwidth h = 0.3 for PG model


k1 = [20 50 100 200 500 1000 2000];

emp =   [37.93  15.42  7.625   3.673   1.543   0.753  0.384];
ex1 =   [31.78  20.48  16.55   13.91   13.15   12.62  12.45];
ex2 =   [45.14  38.99  36.19   34.20  33.98   33.57  33.49];
bino =  [21.59  14.74  11.67   9.842  9.410   9.010  8.914];
poiss = [21.98  17.21  15.63   14.28  14.04   13.74  13.67];
nb =    [34.61  29.08  27.90   26.70  26.30   26.04  25.95];

%

% plot(k1, emp, '-or')
% hold on
% plot(k1, ex1, '-*b')
% plot(k1, ex2, '-*g')
% plot(k1, bino, '-*c')
% plot(k1, poiss, ':k')
% plot(k1, nb, '-.b')


plot(k1, emp, '-or', 'LineWidth', 1.2)
hold on
plot(k1, ex1, '-*b', 'LineWidth', 1.2)
plot(k1, ex2, '-*g', 'LineWidth', 1.2)
plot(k1, bino, '-*c', 'LineWidth', 1.2)
plot(k1, poiss, '--k', 'LineWidth', 1.2)
plot(k1, nb, '-.b', 'LineWidth', 1.2)


legend('Empirical', 'Triangle a=1', 'Triangle a = 2', 'Binomial', 'Poisson', 'Negative Binomial')

% title('Mean of ISE for PG model h = 0.3')






%% bandwidth h = 0.5 for PG model


k1 = [20 50 100 200 500 1000 2000];

emp =   [36.80  14.94  7.687 3.662 1.568 0.753 0.382];
ex1 =   [37.10  27.89  25.90 23.74 22.63 22.59 22.32];
ex2 =   [59.54  54.37  53.61 52.26 51.49 51.68 51.42] ;
bino =  [36.95  29.93  27.80 26.65 25.74 25.95 25.57];
poiss = [26.80  22.84  22.32 21.28 20.75 20.91 20.70];
nb =    [34.81  30.41  29.80 28.57 28.13 28.10 27.96];

%

% plot(k1, emp, '-or')
% hold on
% plot(k1, ex1, '-*b')
% plot(k1, ex2, '-*g')
% plot(k1, bino, '-*c')
% plot(k1, poiss, ':k')
% plot(k1, nb, '-.b')

plot(k1, emp, '-or', 'LineWidth', 1.2)
hold on
plot(k1, ex1, '-*b', 'LineWidth', 1.2)
plot(k1, ex2, '-*g', 'LineWidth', 1.2)
plot(k1, bino, '-*c', 'LineWidth', 1.2)
plot(k1, poiss, '--k', 'LineWidth', 1.2)
plot(k1, nb, '-.b', 'LineWidth', 1.2)


ylim([0 70])

legend('Empirical', 'Triangle a=1', 'Triangle a = 2', 'Binomial', 'Poisson', 'Negative Binomial')

% title('Mean of ISE for PG model h = 0.5')








%% bandwidth h = 0.7 for PG model


k1 = [20 50 100 200 500 1000 2000];

emp =   [37.49 14.94 7.434 3.771 1.503 0.724  0.373   ];
ex1 =   [43.49 34.99 33.59 32.06 31.08 31.08  31.06   ];
ex2 =   [70.99 65.88 65.72 64.73 64.15 64.20  64.29   ] ;
bino =  [65.65 56.81 54.82 53.05 52.53 51.82  52.14   ];
poiss = [36.78 32.93 33.04 32.32 31.92 31.90  32.08   ];
nb =    [38.13 33.89 33.60 32.95 32.43 32.50  32.54   ];

%

% plot(k1, emp, '-or')
% hold on
% plot(k1, ex1, '-*b')
% plot(k1, ex2, '-*g')
% plot(k1, bino, '-*c')
% plot(k1, poiss, ':k')
% plot(k1, nb, '-.b')

plot(k1, emp, '-or', 'LineWidth', 1.2)
hold on
plot(k1, ex1, '-*b', 'LineWidth', 1.2)
plot(k1, ex2, '-*g', 'LineWidth', 1.2)
plot(k1, bino, '-*c', 'LineWidth', 1.2)
plot(k1, poiss, '--k', 'LineWidth', 1.2)
plot(k1, nb, '-.b', 'LineWidth', 1.2)



ylim([0 90])

legend('Empirical', 'Triangle a=1', 'Triangle a = 2', 'Binomial', 'Poisson', 'Negative Binomial')

% title('Mean of ISE for PG model h = 0.7')

















