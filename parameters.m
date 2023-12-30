function parm = parameters

%% All values are annual
parm.gamma = 2; % constant relative risk aversion
parm.alpha = 1/3;  % capital share: alpha; Production function F = K^alpha * L^(1-alpha)
parm.TFP = 1; % total factor productivity is normalized to 1
parm.delta = 0.1;  % Capital depreciation rate: delta
parm.rho = 0.05;   % Time preference parameter(exponential discount rate): rho 
parm.lambda_joblose = 0.4;% job separation intensity, "lambda_1" in the paper
parm.lambda_jobfind_precrisis = 12; % an  initial guess of job finding intensity for computing stationary equilibrium

%% job search cost parameters
% job cost funtion psi(s) = phi*s^(1+kappa)/(1+kappa)
% job finding intensity lambda_2(s) = chi*s
parm.kappa = 0.55;
parm.phi = 0.1; 
parm.chi = 1; % chi is normalized to 1

%% Grid over capital
parm.kmin = 0.01;  % Liquidity constraint
parm.kmax = 30;    % Guess for maximum value of wealth
parm.I = 200;      % number of mesh points on [kmin, kmax]
kmin = parm.kmin;  
kmax = parm.kmax;    
I = parm.I;        
parm.k = linspace(kmin,kmax,I)';  % wealth vector (column vector)
parm.dk = (kmax-kmin)/(I-1);

%% UI benefit parameters
parm.k_low = kmin; % income threshold to receive Ut(unemployment rate) related benefit.
parm.benf_intercept = 0.06*0.3; % "a = 0.018" as in the paper
parm.benf_high = parm.benf_intercept/0.2; % "\bar{y} = 0.09" as in the paper
parm.k_high = parm.k_low+(parm.benf_high-parm.benf_intercept)/0.07; % "k_{H} = 1.039" as in the paper
% acyclical UI benefit policy, as well as recrisis UI befefit policy 
parm.yu_pre = @(k)  parm.benf_intercept+(parm.benf_high - parm.benf_intercept)/(parm.k_high-parm.k_low)*(max(parm.k_low, min(k, parm.k_high))-parm.k_low);

%% parameter in transition dynamics
parm.dt = 0.01; % \Delta_t
parm.T = 25; % terminal time 
parm.M = ceil(parm.T/parm.dt); % number of mesh points on [0, T]
parm.t = linspace(0,parm.T,parm.M)'; % time vector (column vector)
end

