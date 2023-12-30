%unanticipated
function [postCrisis] = dynamics(L_0,parm,preCrisis)

TFP = parm.TFP;

phi = parm.phi;
kappa = parm.kappa;
chi = parm.chi;

alpha = parm.alpha;
delta = parm.delta;
I = parm.I;        
k = parm.k;  %wealth vector (column vector)


dt = parm.dt;
M = parm.M;

lambda_joblose = parm.lambda_joblose;


%Steady state before crisis and at terminal time
v_precrisis = preCrisis.valuefn;
g_precrisis = preCrisis.g;
K_precrisis = preCrisis.K;
L_precrisis = preCrisis.L;
v_T = v_precrisis;
L_T = L_precrisis;
%% Additional job loss at t=0
% L0_0 is employment rate right after crisis;
[g_0] = uniformjobloss(g_precrisis, L_0);

%% Initial guess
K_path=K_precrisis*ones(M,1);
L_path = linspace(L_0,L_T,M)';
w_path = (1-alpha)*TFP.*(K_path.^alpha).*L_path.^(-alpha);
r_path = alpha*TFP.*(K_path.^(alpha-1)).*L_path.^(1-alpha) - delta;

yu_path = zeros(I,M);
tax_path = zeros(M,1);
for i = 1:M
    yu_path(:,i) = parm.benf_calc(k,L_path(i)-L_precrisis);
    tax_path(i) = mean(yu_path(:, i))*(1-L_path(i))/L_path(i);
end
tax_path(1) = g_0(:,1)'*yu_path(:, 1)/L_path(1);
ye_path = w_path - tax_path;


%Preallocation
v_path = zeros(I,2,M);
gg_path = cell(M,1);
g_path = zeros(I,2,M);

Knew_path = zeros(M, 1);
Lnew_path = zeros(M, 1);
r_new_path = zeros(M, 1);
%w_new_path = zeros(M, 1);
yu_new_path = zeros(I, M);
%tax_new_path = zeros(M, 1);

maxit = 500;
convergence_criterion = 10^(-3);
relaxT=0.9;
%%
for it=1:maxit
    fprintf('ITERATION = %d\n',it);
    [Amatrix,v_path,c_path]=HJB_transition(v_T,r_path,ye_path,yu_path,lambda_joblose,parm);
    gg_path{1}=[g_0(:,1); g_0(:,2)];
    for n=1:M-1
        AT=Amatrix{n}';
        %Implicit method in Updating Distribution.
        gg_path{n+1}= (speye(2*I) - AT*dt)\gg_path{n};
        Knew_path(n)=gg_path{n}(1:I)'*k + gg_path{n}(I+1:2*I)'*k;
        Lnew_path(n)= gg_path{n}(I+1:2*I)'*ones(I,1);
        w_path(n) = (1-alpha)*TFP*(Knew_path(n)^alpha)*Lnew_path(n)^(-alpha);
        r_new_path(n) = alpha*TFP*(Knew_path(n)^(alpha-1))*Lnew_path(n)^(1-alpha) - delta;
        yu_new_path(:,n) = parm.benf_calc(k,Lnew_path(n)-L_precrisis);
        tax_path(n) = gg_path{n}(1:I)'*yu_new_path(:, n)/Lnew_path(n);
    end
    Knew_path(M)=gg_path{M}(1:I)'*k + gg_path{M}(I+1:2*I)'*k;
    Lnew_path(M)= gg_path{M}(I+1:2*I)'*ones(I,1);
    w_path(M) = (1-alpha)*TFP*(Knew_path(M)^alpha)*Lnew_path(M)^(-alpha);
    r_new_path(M) = alpha*TFP*(Knew_path(M)^(alpha-1))*Lnew_path(M)^(1-alpha) - delta;
    yu_new_path(:,M) = parm.benf_calc(k,Lnew_path(M)-L_precrisis);
    tax_path(M) = gg_path{M}(1:I)'*yu_new_path(:, M)/Lnew_path(M);
    ye_new_path = w_path - tax_path;
    fprintf('Maximum change in K and L is %.8f \n',max([abs(K_path-Knew_path);abs(L_path-Lnew_path)]));
    if max([abs(K_path-Knew_path);abs(L_path-Lnew_path)])<convergence_criterion
        break
    end
    
    %  Update Kt
    K_path=relaxT.*K_path+(1-relaxT).*Knew_path;
    L_path=relaxT.*L_path+(1-relaxT).*Lnew_path;
    r_path=relaxT.*r_path+(1-relaxT).*r_new_path;
    ye_path=relaxT.*ye_path+(1-relaxT).*ye_new_path;
    yu_path=relaxT.*yu_path+(1-relaxT).*yu_new_path;
end

%%
lbd_jobfind_postcrisis = zeros(I,M);
for i = 1:M
    lbd_jobfind_postcrisis(:,i) = chi.*(chi/phi .*(v_path(:,2,i)-v_path(:,1,i))).^(1/kappa);
end

%%
Ctotal_path = zeros(M, 1);
for i = 1:M
    g_path(:,:,i) = [gg_path{i}(1:I),gg_path{i}(I+1:2*I)];
    Ctotal_path(i) = sum(sum(c_path(:,:,i).*g_path(:,:,i)));
end

postCrisis.valuefn = v_path;
postCrisis.g = g_path;
postCrisis.K_t = K_path;
postCrisis.L_t = L_path;
postCrisis.r_t = r_path;
postCrisis.w_t = w_path;
postCrisis.c_t = c_path;
postCrisis.yu_t = yu_path; % I*M vector
postCrisis.ye_t = ye_path; % M*1 vector
postCrisis.tax_t = tax_path; % NOT tax rate, w_t - ye_t = tax_t = tax rate(tau_t)* w_t
postCrisis.tau_t = postCrisis.tax_t/postCrisis.w_t;
postCrisis.total_expenditure = postCrisis.tax_t.* L_path;

postCrisis.Ctotal = Ctotal_path;
postCrisis.lambda_jobfind = lbd_jobfind_postcrisis;

% figure;
% plot(k, preCrisis.yu,'black', k, postCrisis.yu_t(:,1),'red', k, postCrisis.yu_t(:,M),'green-.')
% title('piecewise linear benefit')
% xlabel('individual income k')
% ylabel('benefit in dollar')
% legend('preCrisis', 't = 0', 't= T')

% figure;
% subplot(3,1,1)
% plot(t, K_path);
% xlabel('Time');
% ylabel('Total capital');
% subplot(3,1,2)
% plot(t, r_path);
% xlabel('Time');
% ylabel('Interest rate');
% subplot(3,1,3)
% plot(t, Ctotal_path);
% xlabel('Time');
% ylabel('Total consumption');

