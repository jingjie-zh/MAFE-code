%% Two different discount rate regimes
parm = parameters;
load("pre.mat")
load("post_ac.mat")
load("post_cc.mat")

load("pre_lesspatient.mat")
load("post_ac_lesspatient.mat")
load("post_cc_lesspatient.mat")

t0 =1;
% low interest case: rho = 0.015
pdf_precrisis_lowinterest = pre_lesspatient.g(:,1)/sum(pre_lesspatient.g(:,1));
cdf_precrisis_lowinterest = cumsum(pdf_precrisis_lowinterest);
L_ac_lowinterest = post_ac_lesspatient.L_t;
U_ac_lowinterest = 1 - L_ac_lowinterest;
lambda_ac_lowinterest = squeeze(post_ac_lesspatient.lambda_jobfind(:,t0));
UDuration_ac_lowinterest = 1./lambda_ac_lowinterest;

L_cc_lowinterest = post_cc_lesspatient.L_t;
U_cc_lowinterest = 1 - L_cc_lowinterest;
lambda_cc_lowinterest = squeeze(post_cc_lesspatient.lambda_jobfind(:,t0));
UDuration_cc_lowinterest = 1./lambda_cc_lowinterest;

excess_urate_lowinterest = U_cc_lowinterest-U_ac_lowinterest;

% high interest case: rho = 0.05
pdf_precrisis_highinterest = preCrisis.g(:,1)/sum(preCrisis.g(:,1));
cdf_precrisis_highinterest = cumsum(pdf_precrisis_highinterest);
L_ac_highinterest = post_ac.L_t;
U_ac_highinterest = 1 - L_ac_highinterest;
lambda_ac_highinterest = squeeze(post_ac.lambda_jobfind(:,t0));
UDuration_ac_highinterest = 1./lambda_ac_highinterest;

L_cc_highinterest = post_cc.L_t;
U_cc_highinterest = 1 - L_cc_highinterest;
lambda_cc_highinterest = squeeze(post_cc.lambda_jobfind(:,t0));
UDuration_cc_highinterest = 1./lambda_cc_highinterest;

excess_urate_highinterest = U_cc_highinterest-U_ac_highinterest;

percentiles = [0.01,0.1,0.25,0.5,0.75,0.9,0.99];
index_highinterest = percentiles;
index_lowinterest = percentiles;
[~,max_i_high] = max(cdf_precrisis_highinterest);
[~,max_i_low] = max(cdf_precrisis_lowinterest);

Xsection_comp = zeros(7,4);
for i = 1:length(percentiles)
    index_highinterest(i) = interp1(cdf_precrisis_highinterest(1:max_i_high),1:max_i_high,percentiles(i),"nearest");
    index_lowinterest(i) = interp1(cdf_precrisis_lowinterest(1:max_i_low),1:max_i_low,percentiles(i),"nearest");
    Xsection_comp(i,:) = [parm.k(index_highinterest(i)),parm.k(index_lowinterest(i)),...
        52*(UDuration_cc_highinterest(index_highinterest(i))-UDuration_ac_highinterest(index_highinterest(i))), ...
        52*(UDuration_cc_lowinterest(index_lowinterest(i))-UDuration_ac_lowinterest(index_lowinterest(i)))];
end
tbl4 = array2table(Xsection_comp, ...
    "VariableNames",["capital_HI","capital_LI","UDuration_HI","UDuration_LI"],...
    'RowNames',{...
    'bottom 1';...
    'bottom 10'; ...
    'bottom 25'; ...
    'median'; ...
    'top 25'; ...
    'top 10'; ...
    'top 1';
    });
display(tbl4)

Tmaxplot = 1;
[max_EU, ind]=max(excess_urate_highinterest-excess_urate_lowinterest);
figure;
plot(parm.t, 100*excess_urate_highinterest, 'b', 'LineWidth', 2); hold on;
plot(parm.t, 100*excess_urate_lowinterest, 'r-.', 'LineWidth', 2); hold on;
xlabel('$t$', 'FontSize',16,'Interpreter','Latex');
ylabel('percent', 'FontSize',16,'Interpreter','Latex');
legend({'$\rho=0.05$', '$\rho=0.015$'},...
    'Location', 'Northeast', 'FontSize',16,'Interpreter','Latex');  legend boxoff
xlim([0 Tmaxplot])
xticks(0:0.25:Tmaxplot)
title('Excess unemployment rate in different economies','FontSize',16,'Interpreter','Latex')
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(h, '-depsc2', 'paperfigures/excessU_comp.eps');

