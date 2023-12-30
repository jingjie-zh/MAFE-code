function [preCrisis]= StationaryE(parm)


%--------Input------------
%TFP: total factor productivity
%lambda: transition rate of employment and unemployment
%at each time, (1,2) vector [lambda1, lambda2]
%--------Output-----------
%v: value function, (I,2) matrix, v(i,j) = v(a_i,y_j)
%g: density function, (I,2)matrix, g(i,j) = g(a_i,y_j)
%K: aggregate capital, which can be computed through g
%L: aggregate labor, which can be computed through g

TFP = parm.TFP;

phi = parm.phi;
kappa = parm.kappa;
chi = parm.chi;

gamma = parm.gamma;
alpha = parm.alpha;
delta = parm.delta;
rho = parm.rho;

kmin = parm.kmin;
kmax = parm.kmax;
I = parm.I;


lambda1 = parm.lambda_jobfind_precrisis;
lambda2 = parm.lambda_joblose;
lambda = [lambda1,lambda2];

k = linspace(kmin,kmax,I)';
dk = (kmax-kmin)/(I-1);
kk = [k,k];


iter_r = 40;
iter_L= 200;
iter_v= 200;
crit_v = 10^(-6);
crit_S = 10^(-5);
crit_L = 10^(-5);
Delta = 1000;

dVf = zeros(I,2);
dVb = zeros(I,2);
c = zeros(I,2);
y = zeros(I,2);
y(:, 1) = parm.yu_pre(k);

rmin = -1*delta;
rmax = 0.999*rho;

S_delta = zeros(1,iter_r);
KS_list = zeros(1,iter_r);
KD_list = zeros(1,iter_r);
r_list = zeros(1,iter_r);
v_delta = zeros(1,iter_v);


% initial guess
r = 0.02;
L = lambda1/(lambda1 + lambda2);
KD = (alpha*TFP/(r + delta))^(1/(1-alpha))*L;
w = (1-alpha)*TFP*KD^alpha*L^(-alpha);
tax = mean(y(:, 1))*(1-L)/L;
y(:,2) = (w-tax)*ones(I,1);

v0(:,1) = (y(:,1) + r.*k).^(1-gamma)/(1-gamma)/rho;
v0(:,2) = (y(:,2) + r.*k).^(1-gamma)/(1-gamma)/rho;
v = v0;

for i=1:iter_r
% update KD and w according r
%     KD = (alpha*TFP/(r + delta))^(1/(1-alpha))*L;
%     w = (1-alpha)*TFP*KD^alpha*L^(-alpha);
% find solutions to HJB with given r, w, denoted by V
% converge when v and V are close enough
    for m = 1:iter_L
%         tau = benifit*(1-L)/L;
%         y = [benifit, 1-tau];
%         yy = ones(I,1)*y;       
        for n=1:iter_v
            V = v;
            % approx of effort and corresponding Phi(effort)
            effort = (chi/phi .*(V(:,2)-V(:,1))).^(1/kappa);
            Phi = phi/(1+kappa)*effort.^(1+kappa);
            % forward difference
            dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/dk;
            dVf(I,:) = (y(I,:) + r.*kmax).^(-gamma); %will never be used, but impose state constraint a<=amax just in case
            % backward difference
            dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/dk;
            dVb(1,:) = (y(1,:) + r.*kmin).^(-gamma); %state constraint boundary condition            
            %consumption and savings with forward difference
            cf = dVf.^(-1/gamma);
            ssf = y + r.*kk - cf;
            %consumption and savings with backward difference
            cb = dVb.^(-1/gamma);
            ssb = y + r.*kk - cb;
            %consumption and derivative of value function at steady state
            c0 = y + r.*kk;
            % dV_upwind makes a choice of forward or backward differences based on
            % the sign of the drift
            If = ssf > 0; %positive drift --> forward difference
            Ib = ssb < 0; %negative drift --> backward difference
            I0 = (1-If-Ib); %at steady state
            c = cf.*If + cb.*Ib + c0.*I0;
            u = c.^(1-gamma)/(1-gamma);            
            %CONSTRUCT MATRIX
            X = -min(ssb,0)/dk;
            Y = -max(ssf,0)/dk + min(ssb,0)/dk;
            Z = max(ssf,0)/dk;
            Aswitch = [-speye(I)*chi.*effort,speye(I)*chi.*effort;speye(I)*lambda(2),-speye(I)*lambda(2)];
            A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
            A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
            A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;          
            if max(abs(sum(A,2)))>10^(-9)
                disp('Improper Transition Matrix')
                break
            end            
            B = (1/Delta + rho)*speye(2*I) - A;
            u_stacked = [u(:,1)-Phi;u(:,2)];
            V_stacked = [V(:,1);V(:,2)];
            b = u_stacked + V_stacked/Delta;
            V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
            V = [V_stacked(1:I),V_stacked(I+1:2*I)];
            Vchange = V - v;
            v = V; % update v
            
            v_delta(n) = max(max(abs(Vchange)));
            if v_delta(n)<crit_v
%                 disp('Value Function Converged, Iteration = ')
%                 disp(n)
                break
            end
            if n == iter_L
                disp('V_fn Not  Converged')
                return
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % FOKKER-PLANCK EQUATION %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        AT = A';
        b = zeros(2*I,1);
        %need to fix one value, otherwise matrix is singular
        i_fix = 1;
        b(i_fix)=.1;
        row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
        AT(i_fix,:) = row;
        %Solve linear system
        gg = AT\b;
        g_sum = gg'*ones(2*I,1);
        gg = gg./g_sum;
        g = [gg(1:I),gg(I+1:2*I)];
        % update L
        DL = g(:,2)'*ones(I,1)-L;
        L = g(:,2)'*ones(I,1);
        % update KD,w,tax,ye according to the latest L and r
        KD = (alpha*TFP/(r + delta))^(1/(1-alpha))*L;
        w = (1-alpha)*TFP*KD^alpha*L^(-alpha);
        tax = g(:,1)'*y(:, 1)/L;
        y(:,2) = (w-tax)*ones(I,1);
        if abs(DL) < crit_L
            disp(['L Converged, Iteration = ',num2str(m)])
            break
        end
        if m == iter_L
                disp('L Not  Converged')
                return
        end
    end
    % update KS according to the latest converged L
    KS = g(:,1)'*k + g(:,2)'*k;
    S_delta(i) = KS - KD;
    KS_list(i) = KS;
    KD_list(i) = KD;
    r_list(i) = r;
    %UPDATE INTEREST RATE
    if S_delta(i)>crit_S
        disp(['Excess Supply, r=',num2str(r), '  S=', num2str(S_delta(i))])
        rmax = r;
        r = 0.5*(r+rmin);
        
    elseif S_delta(i)<-crit_S
        disp(['Excess Demand, r=',num2str(r), '  S=', num2str(S_delta(i))])
        rmin = r;
        r = 0.5*(r+rmax);
    else
        disp('Equilibrium Found')
        break
    end
end
%% Precrisis
K = KS;
preCrisis.Ctotal = sum(c.*g,'all');
preCrisis.lambda_jobfind = chi.*(chi/phi .*(v(:,2)-v(:,1))).^(1/kappa);
preCrisis.valuefn  = v;
preCrisis.K = K;
preCrisis.L = L;
preCrisis.r = alpha* TFP*  K^(alpha-1) * L^(1-alpha);
preCrisis.w = (1-alpha)* TFP *  K^(alpha) * L^(-alpha);
preCrisis.c = c;
preCrisis.g = g;
preCrisis.yu = parm.yu_pre(k); %benefit received by unemploymented (in dollars)
preCrisis.tax = g(:,1)'*preCrisis.yu/L; % NOT tax rate, w - ye = tax = tax rate(tau)* w
preCrisis.ye = preCrisis.w - preCrisis.tax; %wage received by the employmented (in dollars)
preCrisis.tau = preCrisis.tax/preCrisis.w; % tax rate
preCrisis.total_expenditure = preCrisis.tax*L;
preCrisis.kbar_unemployed = g(:,1)'*k/(1-L);
preCrisis.w_2mon = preCrisis.w/6;
%%
% disp(['Interest rate =',num2str(preCrisis.r),', wage = ',num2str(preCrisis.w)])
% disp(['Capital =',num2str(preCrisis.K),', Employment rate = ',num2str(L)])      
% disp(['Total production =',num2str(K^alpha*L^(1-alpha)),', Total consumption = ',num2str(preCrisis.Ctotal)])
%        
