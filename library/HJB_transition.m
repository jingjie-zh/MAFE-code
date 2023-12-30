function [Amatrix_t,v_t,c_t]=HJB_transition(vT,r_t,ye_t,yu_t,lambda_joblose,parm)

%--------Input------------
%vT: boundary condition of value function at T, (I,2) matrix
%K_t: aggregate capital at each time, (N,1) vector
%L_t: aggregate labor at each time, (N,1) vector
%TFP_t: total factor productivity at each time, (N,1) vector
%parm: other parameters

%--------Output-----------
%Amatrix_t: represents the differential operator in HJB equation and it will be
%used in solving KF equation. (N,1) cells, in each cell is a (I*2, I*2)
%matrix
%v_t: value function, (I,2,N) matrix, v(i,j,n) = v(a_i,y_j,t_n)
%c_t: optimal consumption, (I,2,N) matrix, c(i,j,n) = c(a_i,y_j,t_n)

TFP = parm.TFP;

phi = parm.phi;
kappa = parm.kappa;
chi = parm.chi;

gamma = parm.gamma;

rho = parm.rho;
dt = parm.dt;

kmin = parm.kmin;  
kmax = parm.kmax;    
I = parm.I;       

k = linspace(kmin,kmax,I)';
dk = (kmax-kmin)/(I-1);

kk = [k,k];


V = vT;
N=length(r_t);
for n=N:-1:1
    v_t(:,:,n)=V;
    effort = (chi/phi .*(V(:,2)-V(:,1))).^(1/kappa);
    Phi = phi/(1+kappa)*effort.^(1+kappa);
    y = [yu_t(:,n),ones(I,1)*ye_t(n)];
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/dk;
    dVf(I,:) = (y(I,:)  + r_t(n)*kmax).^(-gamma); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/dk;
    dVb(1,:) = (y(1,:)  + r_t(n)*kmin).^(-gamma); %state constraint boundary condition
    %consumption and savings with forward difference
    cf = dVf.^(-1/gamma);
    ssf = y  + r_t(n).*kk - cf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/gamma);
    ssb = y + r_t(n).*kk - cb;
    %consumption and derivative of value function at steady state
    c0 = y + r_t(n).*kk;
    dV0 = c0.^(-gamma);    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift
    If = ssf > 10^(-10); %positive drift --> forward difference
    Ib = ssb < -10^(-10); %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    %make sure backward difference is used at amax
    %Ib(I,:) = 1; If(I,:) = 0;
    %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
    %already taken care of automatically
   
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c = dV_Upwind.^(-1/gamma);
    u = c.^(1-gamma)/(1-gamma);
    c_t(:,:,n)=((dVf+dVb).*(If+Ib)./2 + dV0.*I0).^(-1/gamma);
    c_t(1,:,n)= c(1,:);
    c_t(I,:,n)= c(I,:);
    %CONSTRUCT MATRIX
    X = -min(ssb,0)/dk;
    Y = -max(ssf,0)/dk + min(ssb,0)/dk;
    Z = max(ssf,0)/dk;
    Aswitch = [-speye(I)*chi.*effort,speye(I)*chi.*effort;speye(I)*lambda_joblose,-speye(I)*lambda_joblose];
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    Amatrix_t{n}=A; %save for future reference
    
    B = (1/dt + rho)*speye(2*I) - A;
    u_stacked = [u(:,1)-Phi;u(:,2)];
    V_stacked = [V(:,1);V(:,2)];   
    b = u_stacked + V_stacked/dt;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
end

