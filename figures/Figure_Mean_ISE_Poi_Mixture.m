%%
% the following results are obtained from
% Simulations/PoissonMixture/simulation_ISE

%% bandwidth h = 0.1 for .4*Poi(.5) + .6Poi(10) model


k1 = [20 50 80 100 200 300 400 500 600 700];

emp = [44.02  17.60  11.07  8.986 4.524  3.035 2.192  1.793 1.485   1.309   ];
ex1 = [34.23   14.01 9.039  7.547 3.881  2.776 2.166  1.835 1.600   1.450   ];
bino = [17.41  7.532 5.007  4.202 2.296  1.730 1.413  1.219 1.117   1.036   ];
poiss = [16.88 10.91 9.224  8.530 7.550  7.089 6.873  6.759 6.680   6.655   ];
nb =    [24.35 18.44 16.63  15.83  15.06 14.54  14.30 14.19 14.10   14.12   ];


% plot(k1, emp, '-or')
% hold on
% plot(k1, ex1, '-*b')
% plot(k1, bino, '-*c')
% plot(k1, poiss, ':k')
% plot(k1, nb, '-.b')

plot(k1, emp, '-or', 'LineWidth', 1.2)
hold on
plot(k1, ex1, '-*b', 'LineWidth', 1.2)
plot(k1, bino, '-*c', 'LineWidth', 1.2)
plot(k1, poiss, '--k', 'LineWidth', 1.2)
plot(k1, nb, '-.b', 'LineWidth', 1.2)



legend('Empirical', 'Triangle a=1', 'Binomial', 'Poisson', 'Negative Binomial')

% title('Mean of ISE for .4*Poi(.4) + .6*Poi(10) when h = 0.1')






%% bandwidth h = 0.3 for .4*Poi(.5) + .6Poi(10) model


k1 = [20 50 80 100 200 300 400 500 600 700];

emp = [43.14  16.98  11.28  8.790 4.503  2.943  2.244  1.759  1.509  1.271     ];
ex1 = [26.13  11.61  8.407  7.330 4.770  4.210  3.874  3.524  3.435  3.305     ];
bino = [17.65  8.133  5.816 5.225 3.552  3.199  2.971  2.696  2.728  2.587    ];
poiss = [14.32  9.335 8.446 8.138 7.142  6.974  6.839  6.720  6.715  6.680     ];
nb =    [21.08  16.07 15.46 14.95 13.97  13.70  13.53  13.48  13.39  13.40     ];


% plot(k1, emp, '-or')
% hold on
% plot(k1, ex1, '-*b')
% plot(k1, bino, '-*c')
% plot(k1, poiss, ':k')
% plot(k1, nb, '-.b')


plot(k1, emp, '-or', 'LineWidth', 1.2)
hold on
plot(k1, ex1, '-*b', 'LineWidth', 1.2)
plot(k1, bino, '-*c', 'LineWidth', 1.2)
plot(k1, poiss, '--k', 'LineWidth', 1.2)
plot(k1, nb, '-.b', 'LineWidth', 1.2)


legend('Empirical', 'Triangle a=1', 'Binomial', 'Poisson', 'Negative Binomial')

% title('Mean of ISE for .4*Poi(.4) + .6*Poi(10) when h = 0.3')





%% bandwidth h = 0.5 for .4*Poi(.5) + .6Poi(10) model


k1 = [20 50 80 100 200 300 400 500 600 700];

emp = [45.09  17.73  11.45 9.025  4.538 2.992 2.252  1.771  1.514  1.201  ];
ex1 = [23.448 12.38  9.269 8.589  6.573 5.929 5.836  5.491  5.327  5.159    ];
bino = [22.28  12.38 10.02 9.174  7.572 7.035  6.871 6.651  6.568  6.378     ];
poiss = [13.88 10.60 9.573 9.216  8.646 8.379  8.462 8.316  8.196  8.137     ];
nb =    [19.86 16.46 15.35 14.97  14.32 14.01  14.07 13.92  13.80   13.77    ];


% plot(k1, emp, '-or')
% hold on
% plot(k1, ex1, '-*b')
% plot(k1, bino, '-*c')
% plot(k1, poiss, ':k')
% plot(k1, nb, '-.b')

plot(k1, emp, '-or', 'LineWidth', 1.2)
hold on
plot(k1, ex1, '-*b', 'LineWidth', 1.2)
plot(k1, bino, '-*c', 'LineWidth', 1.2)
plot(k1, poiss, '--k', 'LineWidth', 1.2)
plot(k1, nb, '-.b', 'LineWidth', 1.2)



legend('Empirical', 'Triangle a=1', 'Binomial', 'Poisson', 'Negative Binomial')

% title('Mean of ISE for .4*Poi(.4) + .6*Poi(10) when h = 0.5')






%% bandwidth h = 0.7 for .4*Poi(.5) + .6Poi(10) model


k1 = [20 50 80 100 200 300 400 500 600 700];

emp = [43.74  17.57   11.26   8.954  4.407 3.020  2.234  1.751  1.540  1.295        ];
ex1 = [22.62   13.16  10.54   9.676  8.134 7.698  7.502  7.064  7.241  7.068        ];
bino = [32.78   19.76  17.10  15.63  14.05 13.52  13.11  12.72  12.99  12.64        ];
poiss = [15.68  12.75  11.87  11.69  11.26 11.22  11.18  10.91  11.14   11.00      ];
nb =    [19.69  16.96  15.96  15.91  15.27 15.20  15.15  14.94   15.05  14.99       ];


% plot(k1, emp, '-or')
% hold on
% plot(k1, ex1, '-*b')
% plot(k1, bino, '-*c')
% plot(k1, poiss, ':k')
% plot(k1, nb, '-.b')

plot(k1, emp, '-or', 'LineWidth', 1.2)
hold on
plot(k1, ex1, '-*b', 'LineWidth', 1.2)
plot(k1, bino, '-*c', 'LineWidth', 1.2)
plot(k1, poiss, '--k', 'LineWidth', 1.2)
plot(k1, nb, '-.b', 'LineWidth', 1.2)



legend('Empirical', 'Triangle a=1', 'Binomial', 'Poisson', 'Negative Binomial')

% title('Mean of ISE for .4*Poi(.4) + .6*Poi(10) when h = 0.7')


