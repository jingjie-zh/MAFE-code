
%% Plots: All as a function of initial shock size u0
parm = parameters;
load("CompShockSize.mat")

[~, tn] = min(abs(parm.t-1));
t_index = 18;

UIoverGDP_cc = sum(UIexpense_cc_path(1:tn,:),1)./sum(Y_cc_path(1:tn,:),1);
UIoverGDP_ac = sum(UIexpense_ac_path(1:tn,:),1)./sum(Y_ac_path(1:tn,:),1);
excess_urate = L_ac_path(t_index,:) - L_cc_path(t_index,:);

figure;

subplot(1,3,1)
plot(100*u0_grid, 100*UIoverGDP_cc, 'b', 'LineWidth', 2); hold on;
plot(100*u0_grid, 100*UIoverGDP_ac, 'k--', 'LineWidth', 2); hold on;
xlabel('$U_0$', 'Interpreter', 'Latex', 'FontSize',14);
legend({'$\eta = 10$','$\eta = 0$'},...
       'Location', 'Northwest', 'FontSize',14,'Interpreter','Latex');  
legend boxoff
axis('square');
xlim([6 10])
xticks(6:2:10)
ylim([0.3 1.1])
yticks(0.3:0.2:1.1)
box off
title('A. Total UI/GDP in a year', 'Interpreter', 'Latex', 'FontSize',12);

subplot(1,3,2)
slope = (excess_urate(end-1)-excess_urate(end))/(u0_grid(end-1)-u0_grid(end));
linear_comp = 100*(excess_urate(end)+slope*(u0_grid-0.06));
plot(100*u0_grid, 100*excess_urate, 'b', 'LineWidth', 2); hold on;
plot(100*u0_grid, linear_comp, 'k--', 'LineWidth', 2); hold on;
xlabel('$U_0$', 'Interpreter', 'Latex', 'FontSize',14);
ylabel('percent', 'Interpreter', 'Latex', 'FontSize',14);
axis('square');
xlim([6 10])
xticks(6:2:10)
ylim([0 1])
yticks(0:0.2:1)
box off
title('B. Excess unemployment, $t=2$ months', 'Interpreter', 'Latex', 'FontSize',12);

subplot(1,3,3)
plot(100*u0_grid, 100*(Y_cc_path(t_index,:)-Y_precrisis)/Y_precrisis, 'b', 'LineWidth', 2); hold on;
plot(100*u0_grid, 100*(Y_ac_path(t_index,:)-Y_precrisis)/Y_precrisis, 'k--', 'LineWidth', 2); 
legend({'$\eta = 10$','$\eta = 0$'},...
       'Location', 'Southwest', 'FontSize',14,'Interpreter','Latex');  
legend boxoff
xlabel('$U_0$', 'Interpreter', 'Latex', 'FontSize',14);
ylabel('percent', 'Interpreter', 'Latex', 'FontSize',14);
axis('square');
xlim([6 10])
xticks(6:2:10)
ylim([-1.8 0])
yticks(-1.8:0.2:0)
box off
title('C. Output, $t=2$ months', 'Interpreter', 'Latex', 'FontSize',12);


h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(h, '-depsc2', 'paperfigures/compstats.eps');

