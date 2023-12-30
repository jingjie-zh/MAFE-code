%% Aggregate results
parm = parameters;
load("pre.mat")
load("post_ac.mat")
load("post_cc.mat")

L_ac_U = post_ac.L_t;
Y_ac_U = parm.TFP.*(post_ac.K_t.^parm.alpha).*(post_ac.L_t).^(1-parm.alpha);
w_ac_U = post_ac.w_t;

L_cc_U = post_cc.L_t;
Y_cc_U = parm.TFP.*(post_cc.K_t.^parm.alpha).*(post_cc.L_t).^(1-parm.alpha);
w_cc_U = post_cc.w_t;

Tmaxplot = 1;

figure;
subplot(1,3,1)
plot(parm.t, 100*(1-L_ac_U), 'b', 'LineWidth', 2); hold on;
plot(parm.t, 100*(1-L_cc_U), 'r-.', 'LineWidth', 2);
xlabel('$t$', 'FontSize',12,'Interpreter','Latex');
ylabel('percent', 'FontSize',12,'Interpreter','Latex');
legend({'$\eta=0$', '$\eta=10$'},...
    'Location', 'Northeast', 'FontSize',12,'Interpreter','Latex'); legend boxoff
title('A: Unemployment rate','FontSize',12,'Interpreter','Latex', 'FontSize',12)
xlim([0 Tmaxplot])
xticks(0:0.25:Tmaxplot)
ylim([5 12])
yticks(5:2:12)
box off
axis square

subplot(1,3,2)
plot(parm.t, 100*(Y_ac_U - Y_ac_U(1))/Y_ac_U(1), 'b', 'LineWidth', 2); hold on;
plot(parm.t, 100*(Y_cc_U - Y_cc_U(1))/Y_cc_U(1), 'r-.', 'LineWidth', 2);
xlabel('$t$', 'FontSize',12,'Interpreter','Latex');
ylabel('percent', 'FontSize',12,'Interpreter','Latex');
legend({'$\eta=0$', '$\eta=10$'},...
    'Location', 'Southeast', 'FontSize',12,'Interpreter','Latex');  legend boxoff
title('B: Output growth','FontSize',12,'Interpreter','Latex', 'FontSize',12)
xlim([0 Tmaxplot])
xticks(0:0.25:Tmaxplot)
ylim([0 4])
yticks(0:1:4)
box off
axis square

subplot(1,3,3)
plot(parm.t, 100*(w_ac_U - w_ac_U(1))/w_ac_U(1), 'b', 'LineWidth', 2); hold on
plot(parm.t, 100*(w_cc_U - w_cc_U(1))/w_cc_U(1), 'r-.', 'LineWidth', 2); hold on

xlabel('$t$', 'FontSize',12,'Interpreter','Latex');
ylabel('percent', 'FontSize',12,'Interpreter','Latex');
legend({'$\eta=0$', '$\eta=10$'},...
    'Location', 'Northeast', 'FontSize',12,'Interpreter','Latex');  legend boxoff
xlim([0 Tmaxplot])
xticks(0:0.25:Tmaxplot)
ylim([-2 0.5])
yticks(-2: 0.5: 0.5)
title('C: Wage change','FontSize',12,'Interpreter','Latex', 'FontSize',12)
box off
axis square

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(h, '-depsc2', 'paperfigures/aggregateresults_baseline.eps');


