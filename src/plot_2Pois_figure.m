function plot_2Pois_figure(param)
    % this is the plot for two component Poisson mixture 0.4*Poi(0.5) +
    % 0.6*Poi(10), the table title is Estimates Using Different Kernels for $0.4 Poi(0.5) + 0.6Poi(10)$
    % plot_table_onefigure('pi1' | 'lambda1' | 'lambda2')
    % Default: 'pi1'
    if nargin==0, param = 'pi1'; end
    param = lower(param);
    assert(ismember(param, {'pi1','lambda1','lambda2'}), ...
        'param must be ''pi1'', ''lambda1'', or ''lambda2''.');

    % Sample sizes
    n = [20 50 100 200];

    methods = {'EM','HD','VNED'};
    kernels = {'Empirical','Triangle_a1','Triangle_a2','Binomial','Poisson','NegBinomial'};

    % ---------- TRUE VALUES ----------
    truths.pi1     = 0.4;
    truths.lambda1 = 0.5;
    truths.lambda2 = 10;

    % ---------- DATA: averages only (columns: EM, HD, VNED) ----------
    % Empirical
    data.Empirical.pi1     = [0.409 0.507 0.449; 0.397 0.435 0.423; 0.399 0.417 0.415; 0.398 0.407 0.408];
    data.Empirical.lambda1 = [0.544 0.410 0.456; 0.490 0.442 0.463; 0.502 0.478 0.482; 0.502 0.489 0.488];
    data.Empirical.lambda2 = [9.904 9.605 9.556; 10.01 9.788 9.866; 10.02 9.849 9.917; 10.07 9.879 9.941];

    % Triangle a=1
    data.Triangle_a1.pi1     = [0.409 0.488 0.413; 0.397 0.406 0.390; 0.399 0.398 0.395; 0.398 0.397 0.398];
    data.Triangle_a1.lambda1 = [0.544 0.487 0.589; 0.490 0.570 0.621; 0.502 0.584 0.600; 0.502 0.551 0.553];
    data.Triangle_a1.lambda2 = [9.904 9.604 9.639; 10.01 9.772 9.873; 10.02 9.841 9.922; 10.07 9.875 9.942];

    % Triangle a=2
    data.Triangle_a2.pi1     = [0.409 0.467 0.391; 0.397 0.401 0.385; 0.399 0.397 0.394; 0.398 0.396 0.397];
    data.Triangle_a2.lambda1 = [0.544 0.498 0.631; 0.490 0.572 0.620; 0.502 0.584 0.598; 0.502 0.553 0.554];
    data.Triangle_a2.lambda2 = [9.904 9.578 9.638; 10.01 9.768 9.873; 10.02 9.835 9.915; 10.07 9.869 9.934];

    % Binomial
    data.Binomial.pi1     = [0.409 0.517 0.425; 0.397 0.389 0.368; 0.399 0.363 0.357; 0.398 0.350 0.350];
    data.Binomial.lambda1 = [0.544 0.354 0.408; 0.490 0.438 0.463; 0.502 0.492 0.502; 0.502 0.513 0.515];
    data.Binomial.lambda2 = [9.904 9.441 9.453; 10.01 9.503 9.554; 10.02 9.510 9.562; 10.07 9.491 9.547];

    % Poisson
    data.Poisson.pi1     = [0.409 0.565 0.473; 0.397 0.466 0.447; 0.399 0.451 0.448; 0.398 0.448 0.449];
    data.Poisson.lambda1 = [0.544 0.426 0.519; 0.490 0.604 0.651; 0.502 0.765 0.777; 0.502 0.878 0.861];
    data.Poisson.lambda2 = [9.904 9.380 9.396; 10.01 9.582 9.649; 10.02 9.752 9.818; 10.07 9.911 9.942];

    % Negative Binomial
    data.NegBinomial.pi1     = [0.409 0.636 0.563; 0.397 0.532 0.521; 0.399 0.517 0.519; 0.398 0.517 0.520];
    data.NegBinomial.lambda1 = [0.544 0.466 0.565; 0.490 0.681 0.732; 0.502 0.903 0.912; 0.502 1.089 1.055];
    data.NegBinomial.lambda2 = [9.904 9.400 9.366; 10.01 9.586 9.618; 10.02 9.812 9.832; 10.07 10.09 10.05];

    % ---------- Plot (3x2 grid using subplot) ----------
    figure('Color','w');
    colors = lines(3); markers = {'o','s','^'};

    for k = 1:numel(kernels)
        subplot(3,2,k); hold on;
        M = data.(kernels{k}).(param);  % [4 x 3], columns: EM, HD, VNED
        for j = 1:3
            plot(n, M(:,j), '-', 'Color', colors(j,:), ...
                'LineWidth', 1.8, 'Marker', markers{j}, 'MarkerSize', 6);
        end

        % Add horizontal line for truth
        yline(truths.(param), 'k--', 'LineWidth', 2);

        title(strrep(kernels{k},'_',' '));
        grid on; box on;
        if k > 4, xlabel('Sample size n'); end
        if mod(k,2)==1
            switch param
                case 'pi1',     ylabel('Average \pi_1');
                case 'lambda1', ylabel('Average \lambda_1');
                case 'lambda2', ylabel('Average \lambda_2');
            end
        end
    end

    % One legend for all subplots
    hL = legend(methods, 'Orientation', 'horizontal');
    set(hL, 'Position', [0.35 0.01 0.3 0.05]);  % adjust if needed

    % Overall title
    switch param
        case 'pi1',     ttl = 'Average Estimates of \pi_1 Across Kernels';
        case 'lambda1', ttl = 'Average Estimates of \lambda_1 Across Kernels';
        case 'lambda2', ttl = 'Average Estimates of \lambda_2 Across Kernels';
    end
    sgtitle(ttl);  % available in R2019a
end
