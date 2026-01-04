% written on 06/02/2019

% this is to plot RAF for mixture model paper




%%

syms x
fplot(x ,[-1 6], '-r')
hold on
fplot(2.*((x+1).^.5-1 ), [-1 6], '--b')
fplot(2-2.*exp(-x)-x.*exp(-x), [-1 6], ':k')
fplot(exp(-1./(1+x)+1) -1, [-1 6], '-.g')

legend('LD', 'HD', 'NED', 'VNED' )

xlabel('\delta')
ylabel('A(\delta)')

% title('Residual Adjustment Function A(\delta)')

%%
syms x
clf

% soft (pastel-ish) colors
c.red   = [0.80 0.40 0.40];
c.blue  = [0.40 0.60 0.85];
c.green = [0.45 0.75 0.60];
c.gray  = [0.45 0.45 0.45];

h1 = fplot(x, [-1 6], 'LineWidth', 2.5, 'Color', c.red); hold on
h2 = fplot(2*((x+1).^(1/2) - 1), [-1 6], 'LineStyle','--', 'LineWidth', 2.5, 'Color', c.blue);
h3 = fplot(2 - 2*exp(-x) - x.*exp(-x), [-1 6], 'LineStyle',':',  'LineWidth', 2.5, 'Color', c.gray);
h4 = fplot(exp(-1./(1+x) + 1) - 1,     [-1 6], 'LineStyle','-.', 'LineWidth', 2.5, 'Color', c.green);

legend([h1 h2 h3 h4], 'LD', 'HD', 'NED', 'VNED', 'Location','northwest'); % top-left
xlabel('\delta'); ylabel('A(\delta)');

grid on; box on;
set(gca,'LineWidth',1.2)   % thicker axes

