%% Cross-sectional effects
parm = parameters;
load("pre.mat")
load("post_ac.mat")
load("post_cc.mat")


t0 = 1;
c_preCrisis_E = preCrisis.c(:,2);
lambda_ac_U = squeeze(post_ac.lambda_jobfind(:,t0));
UD_ac = 1./lambda_ac_U;
c_ac_U = squeeze(post_ac.c_t(:,1,t0));
FCC_ac = (c_ac_U-c_preCrisis_E)./c_preCrisis_E;
VF_ac = squeeze(post_ac.valuefn(:,1,t0));

lambda_cc_U = squeeze(post_cc.lambda_jobfind(:,t0));
UD_cc = 1./lambda_cc_U;
c_cc_U = squeeze(post_cc.c_t(:,1,t0));
FCC_cc = (c_cc_U-c_preCrisis_E)./c_preCrisis_E;
VF_cc = squeeze(post_cc.valuefn(:,1,t0));

g = sum(preCrisis.g,2);
G = g;
for i = 2:length(g)
    G(i) = g(i)+G(i-1);
end
[~,max_i] = max(G);
% find out index of different percentile
i_top01 = interp1(G(1:max_i),1:max_i,0.99,"nearest");
i_top10 = interp1(G(1:max_i),1:max_i,0.9,"nearest");
i_top25 = interp1(G(1:max_i),1:max_i,0.75,"nearest");
i_bot50 = interp1(G(1:max_i),1:max_i,0.5,"nearest");
i_bot25 = interp1(G(1:max_i),1:max_i,0.25,"nearest");
i_bot10 = interp1(G(1:max_i),1:max_i,0.10,"nearest");
i_bot05 = interp1(G(1:max_i),1:max_i,0.05,"nearest");
i_bot01 = interp1(G(1:max_i),1:max_i,0.01,"nearest");

Xsection = zeros(7,6);
indexList = [i_bot01,i_bot10,i_bot25,i_bot50, i_top25, i_top10,i_top01];
for n = 1:7
    i = indexList(n);
    Xsection(n,:) = [52*UD_ac(i), 52*UD_cc(i), ...
        100*FCC_ac(i), 100*FCC_cc(i),...
        VF_ac(i), VF_cc(i)];
end

tbl3 = array2table(Xsection, ...
    "VariableNames",["UD_ac","UD_cc","FCC_ac","FCC_cc","VF_ac","VF_cc"],...
    'RowNames',{...
    'bottom 1';...
    'bottom 10'; ...
    'bottom 25'; ...
    'median'; ...
    'top 25'; ...
    'top 10'; ...
    'top 1';
    });
display(tbl3)