% This is Lema Image Analysis with Poisson Mixture

Y_c = rgb2gray(image_c);
Y = double(Y_c);

Y = Y(:);

%% EM 2 component
[EMIC2, hFig, params, pi_hat, lambda_hat] = em_bic_plot( ...
    Y, 2, 'LabelLevels',[1 150], 'AutoImage', true, ...
    'MaxIter', 800, 'Tol', 1e-8, 'Init', 'quantile');

EMIC2.BIC

%% EM 3 component
[EMIC3, hFig, params, pi_hat, lambda_hat] = em_bic_plot( ...
     Y, 3, 'LabelLevels',[1 100 200], 'AutoImage', true, ...
    'MaxIter', 800, 'Tol', 1e-8, 'Init', 'quantile');

EMIC3.BIC

%% EM 4 component
[EMIC4, hFig, params, pi_hat, lambda_hat] = em_bic_plot( ...
     Y, 4, 'LabelLevels',[1 100 150 225], 'AutoImage', true, ...
    'MaxIter', 800, 'Tol', 1e-8, 'Init', 'quantile');

EMIC4.BIC

%% EM 5 component
[EMIC5, hFig, params, pi_hat, lambda_hat] = em_bic_plot( ...
    Y, 5, 'LabelLevels', [1 50 150 175 225], 'AutoImage', true, ...
    'MaxIter', 800, 'Tol', 1e-8, 'Init', 'quantile');

EMIC5.BIC

%% HMIX 2 component
[HIC2, hFig, params, pi_hat, lambda_hat] = hmix_hic_plot( ...
    Y, 2, 'LabelLevels',[1 150], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);


%% HMIX 3 component
[HIC3, hFig, params, pi_hat, lambda_hat] = hmix_hic_plot( ...
    Y, 3, 'LabelLevels',[1 100 200], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);



%% HMIX 4 component
[HIC4, hFig, params, pi_hat, lambda_hat] = hmix_hic_plot( ...
    Y, 4, 'LabelLevels',[1 100 150 225], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);


%% HMIX 5 component

[HIC5, hFig, params, pi_hat, lambda_hat] = hmix_hic_plot( ...
    Y, 5, 'LabelLevels',[1 50 150 175 225], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);



%% VNEDMIX 2 component
[VNEDIC2, hFig, params, pi_hat, lambda_hat] = vnedmix_dic_plot( ...
    Y, 2, 'LabelLevels',[1 150], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);


%% VNEDMIX 3 component
[VNEDIC3, hFig, params, pi_hat, lambda_hat] = vnedmix_dic_plot( ...
    Y, 3, 'LabelLevels',[1 100 200], 'AutoImage', true, ...
    'Init','quantile', 'MaxIter',800, 'Tol',1e-8);


%% VNEDMIX 4 component
[VNEDIC4, hFig, params, pi_hat, lambda_hat] = vnedmix_dic_plot( ...
    Y, 4, 'LabelLevels',[1 100 150 225], 'AutoImage', true, ...
    'MaxIter',800, 'Tol',1e-8, 'Init','quantile');

%% VNEDMIX 5 component
[VNEDIC5, hFig, params, pi_hat, lambda_hat] = vnedmix_dic_plot( ...
    Y, 5, 'LabelLevels',[1 50 150 175 225], 'AutoImage', true, ...
    'MaxIter',800, 'Tol',1e-8, 'Init','quantile');
