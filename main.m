clear all;
close all;
clc;
addpath library
% pre-crisis stationary equilibrium
parm = parameters;
[preCrisis]= StationaryE(parm);
% post-crisis time-dependent equilibrium 
L_0 = 0.9;
% countercyclical policy
parm.eta = -10;
parm.benf_calc = @(k,Delta_L) parm.eta*Delta_L + parm.yu_pre(k);
[post_cc] = dynamics(L_0,parm,preCrisis);
% acyclical policy
parm.eta = 0;
parm.benf_calc = @(k,Delta_L) parm.eta*Delta_L + parm.yu_pre(k);
[post_ac] = dynamics(L_0,parm,preCrisis);

save("pre.mat","preCrisis")
save("post_ac.mat","post_ac")
save("post_cc.mat","post_cc")
%%
% reset rho from 0.05 to 0.015
parm = parameters;
parm.rho = 0.015;
% pre-crisis stationary equilibrium
[pre_lesspatient]= StationaryE(parm); 
% post-crisis time-dependent equilibrium 
L_0 = 0.9;
% countercyclical policy
parm.eta = -10;
parm.benf_calc = @(k,Delta_L) parm.eta*Delta_L + parm.yu_pre(k);
[post_cc_lesspatient] = dynamics(L_0,parm,pre_lesspatient);
% acyclical policy
parm.eta = 0;
parm.benf_calc = @(k,Delta_L) parm.eta*Delta_L + parm.yu_pre(k);
[post_ac_lesspatient] = dynamics(L_0,parm,pre_lesspatient);

save("pre_lesspatient.mat","pre_lesspatient")
save("post_ac_lesspatient.mat","post_ac_lesspatient")
save("post_cc_lesspatient.mat","post_cc_lesspatient")
%% comparativestatics: excess unemployment vs shock size
% reset parameters
parm = parameters;
% same stationary equilbirum for all cases
[preCrisis]= StationaryE(parm);
Y_precrisis = parm.TFP*((preCrisis.K)^parm.alpha).*(preCrisis.L)^(1-parm.alpha);
u_precrisis = 1 - preCrisis.L;

gridN = 9;
Lrategrid = linspace(0.9, 0.94, gridN);
u0_grid = 1 - Lrategrid;
% Preallocation
UIexpense_cc_path = nan(parm.M, gridN);
UIexpense_ac_path = nan(parm.M, gridN);
Y_cc_path = nan(parm.M, gridN);
Y_ac_path = nan(parm.M, gridN);
L_cc_path = nan(parm.M, gridN);
L_ac_path = nan(parm.M, gridN);
K_cc_path = nan(parm.M, gridN);
K_ac_path = nan(parm.M, gridN);

%%%%%%%%%% countercyclical 
parm.eta = -10;
eta_countercyclical = parm.eta;
parm.benf_calc = @(k,Delta_L) parm.eta*Delta_L + parm.yu_pre(k);

for i = 1:gridN
    Lnow = Lrategrid(i);
    [postCrisis] = dynamics(Lnow, parm, preCrisis); 
    UIexpense_cc_path(:,i) = postCrisis.total_expenditure;
    %reminder: postCrisis.total_expenditure = postCrisis.tax_t.* L_path 
    % = sum(postCrisis.g(:,1,t0).*postCrisis{i}.yu_t(:,t0));
    L_cc_path(:,i)= postCrisis.L_t;
    K_cc_path(:,i)= postCrisis.K_t;
    Y_cc_path(:,i)=parm.TFP.*((postCrisis.K_t).^parm.alpha).*(postCrisis.L_t).^(1-parm.alpha);
end
%%%%%%%% acyclical 
parm.eta = 0;
parm.benf_calc = @(k,Delta_L) parm.eta*Delta_L + parm.yu_pre(k);
eta_acyclical = parm.eta;

for i = 1:gridN
    Lnow = Lrategrid(i);
    [postCrisis] = dynamics(Lnow, parm, preCrisis); 
    UIexpense_ac_path(:,i) = postCrisis.total_expenditure;%
    L_ac_path(:,i)= postCrisis.L_t;
    K_ac_path(:,i)= postCrisis.K_t;
    Y_ac_path(:,i)=parm.TFP.*((postCrisis.K_t).^parm.alpha).*(postCrisis.L_t).^(1-parm.alpha);
end

save("CompShockSize.mat",...
    "UIexpense_cc_path",...
    "UIexpense_ac_path",...
    "Y_cc_path",...
    "Y_ac_path",...
    "L_cc_path",...
    "L_ac_path",...
    "L_cc_path",...
    "L_ac_path",...
    "u0_grid",...
    "Y_precrisis")

