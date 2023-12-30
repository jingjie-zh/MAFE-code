%% Shock and U.I. policy
parm = parameters;
load("pre.mat")

gprecrisis = preCrisis.g;
L0 = 0.9;
gpostcrisis = uniformjobloss(gprecrisis, L0);


figure;
subplot(1,3,1)
eps_i = 1;
plot(parm.k, (squeeze(gprecrisis(:,eps_i))), 'k', 'LineWidth', 2); hold on;
plot(parm.k, (squeeze(gpostcrisis(:,eps_i))), 'r-.', 'LineWidth', 2);
xlabel('$k$', 'Interpreter', 'Latex', 'FontSize',14);
ylabel('$g_1$', 'Interpreter', 'Latex', 'FontSize',14);
axis('square');
xlim([0 10])
xticks(0:2:10)
ylim([0 0.006])
yticks(0:0.002:0.006);
box off
legend({'$t<0$','t=0'},...
    'Location', 'Northeast', 'FontSize',14,'Interpreter','Latex');  legend boxoff
legend boxoff
title('A. Unemployed', 'Interpreter', 'Latex', 'FontSize',12);

subplot(1,3,2)
eps_i = 2;
plot(parm.k, (squeeze(gprecrisis(:,eps_i))), 'k', 'LineWidth', 2); hold on;
plot(parm.k, (squeeze(gpostcrisis(:,eps_i))), 'r-.', 'LineWidth', 2);

xlabel('$k$', 'Interpreter', 'Latex', 'FontSize',14);
ylabel('$g_2$', 'Interpreter', 'Latex', 'FontSize',14);
axis('square');
xlim([0 10])
xticks(0:2:10)
ylim([0 0.06])
yticks(0:0.02:0.06)
box off
legend({'$t<0$','t=0'},...
    'Location', 'Northeast', 'FontSize',14,'Interpreter','Latex');  legend boxoff
legend boxoff
title('B. Employed', 'Interpreter', 'Latex', 'FontSize',12);

% U.I. policies
parm.eta = -10;
parm.benf_calc = @(k,Delta_L) parm.eta*Delta_L + parm.yu_pre(k);
postcrisis_yu = parm.benf_calc(parm.k,L0-preCrisis.L);

subplot(1,3,3)
plot(parm.k, preCrisis.yu, 'k', 'LineWidth', 2); hold on;
plot(parm.k, postcrisis_yu, 'r-.', 'LineWidth', 2);
xlabel('$k$', 'Interpreter', 'Latex', 'FontSize',14);
ylabel('$y_1$', 'Interpreter', 'Latex', 'FontSize',14);
axis('square');
xlim([0 4])
xticks(0:1:4)
ylim([0 0.8])
yticks(0:0.2:0.8);
box off
legend({'$t<0$','t=0'},...
    'Location', 'Northeast', 'FontSize',14,'Interpreter','Latex');  legend boxoff
legend boxoff
title('C. UI benefit', 'Interpreter', 'Latex', 'FontSize',12);


h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(h, '-depsc2', 'paperfigures/shock_baseline.eps');
