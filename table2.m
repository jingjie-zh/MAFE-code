%% calibration targets
parm = parameters;
load("pre.mat")
% average elasticity of unemployment duration with respect to unemployment benefits
ind = find(parm.k>parm.k_high,1)-1;
part1 =diff(1./preCrisis.lambda_jobfind(1:ind))./diff(preCrisis.yu(1:ind));
part2=preCrisis.yu(1:ind-1).*preCrisis.lambda_jobfind(1:ind-1);
average_Uduration_elasticity =  sum(part1.*part2.*preCrisis.g(1:ind-1,1))/sum(preCrisis.g(1:ind-1,1));

% incomesensitivity_UI
totalincome = preCrisis.r*parm.k + preCrisis.w;
k2 = sum(parm.k < parm.k_high);
incomesensitivity_UI = (parm.yu_pre(parm.k(k2)) - parm.yu_pre(parm.k(1)))/(totalincome(k2) - totalincome(1));

u_pdf_precrisis = preCrisis.g(:,1)/sum(preCrisis.g(:,1));
meanincome = sum(totalincome.*u_pdf_precrisis);
threshold_income = preCrisis.r*parm.k_high + preCrisis.w;
minUI_over_maxUI = parm.benf_intercept/parm.benf_high; % Data value from Krueger Meyer (2002) review


% Mean and Median replacement rates from Ganong, Noel, Vavra (J. of Public
% Economics, 2020)
% incomesensitivity_UI from Ganong, Noel, Vavra (J. of Public Economics,
% 2020) estimate for Nevada which is in the middle of the US distribution.
%Elasticity of unemployment duration with respect to the benefit level is from Chetty

model = round([incomesensitivity_UI; meanincome; threshold_income; minUI_over_maxUI; average_Uduration_elasticity;],4);
data = [0.52; 886; 902; 0.2; 0.54];

tbl0 = table( model,  data, 'RowNames',...
    {...
    'UI benefit sensitivity to total income';...
    'Mean total income';...
    'Income at kink';...
    'Min UI divided by max UI';...
    '(average) Elasticity of U. duration w.r.t. b';...
    });
display(tbl0)
%% Table: Precrisis values

Y0 = parm.TFP*(preCrisis.K^parm.alpha)*(preCrisis.L)^(1-parm.alpha);
aggregate = [ preCrisis.L; preCrisis.w; preCrisis.K; preCrisis.r; preCrisis.Ctotal; Y0; preCrisis.total_expenditure/Y0];
tbl2 = table( aggregate, 'RowNames',...
    {...
    'L'
    'wage'; ...
    'K'; ...
    'rental rate'; ...
    'Aggregate C'; ...
    'Aggregate Output'; ...
    'Total UI/GDP';
    });
display(tbl2)

