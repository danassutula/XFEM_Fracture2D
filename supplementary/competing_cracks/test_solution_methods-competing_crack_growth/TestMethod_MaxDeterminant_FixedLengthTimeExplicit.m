home;

% Goal: test if 'MethodOfMaximumDeterminat' is adequate for determining the
% critical crack tips assuming a uniform extension of the critical tips.
%
% The algorithm is robust in the following cases:
%   1. non-competing crack growth (Gs_i ~= Gs_j), negative-definite, 
% positive-definite or indefinite Hs.
%   2. competing crack growth (Gs_i == Gs_j), negative-definite or
% positive-definite Hs (but not for indefinite Hs).
%
% The algorithm will not be generally robust for the case of:
%   1. competing crack growth (Gs_i == Gs_j) and an indefinite Hs.

if ~exist('flag_loadCase','var') 
    flag_loadCase = 0;
end

while 1

%% INPUTS
if flag_loadCase == 0
    
    %%% CHOOSE:
    da_tot = 1.0 % total fracture growth
    n_tips = 12 % n. crack tips participating
    n_mesh = 10 % n_disc-1 many times to half the crack extension length da_inc
    
    %%% CHOOSE:
    type_Gs = 'uniform';
%     type_Gs = 'nonunif';
    
    %%% CHOOSE:
%     type_Hs = 'indefn';  % is indefinite
%     type_Hs = 'posdef'; % is positive def.
    type_Hs = 'negdef'; % is negative def.
    
    NAME = sprintf('GsType(%s)_HsType(%s)_nTips(%i)',...
        type_Gs,type_Hs,n_tips);
    
end


%% LOAD/GENERATE A CASE 
if flag_loadCase == 1
    
    n_tips = INPUTS.n_tips;
    n_mesh = INPUTS.n_mesh;
    da_tot = INPUTS.da_tot;
    
    type_Gs = INPUTS.type_Gs;
    type_Hs = INPUTS.type_Hs;

    Gs = INPUTS.Gs;
    Hs = INPUTS.Hs;
    
elseif flag_loadCase == 0
    
    INPUTS = [];

    INPUTS.da_tot = da_tot;
    INPUTS.n_tips = n_tips;
    INPUTS.n_mesh = n_mesh;
    
    INPUTS.type_Gs = type_Gs;
    INPUTS.type_Hs = type_Hs;
    
    Gs = rand(n_tips,1);
    
    switch type_Gs
        case 'uniform'; Gs(:) = Gs(1);
        case 'nonunif'; % do nothing
    end
    
    % random matrix from
    Hs = rand(n_tips,n_tips);
    
    switch type_Hs
        case 'indefn'; Hs = Hs'*diag(2*(rand(n_tips,1)-0.5))*Hs;
        case 'posdef'; Hs = Hs'*diag(rand(n_tips,1))*Hs;
        case 'negdef'; Hs = Hs'*diag(-rand(n_tips,1))*Hs;
    end
    
    INPUTS.Gs = Gs;
    INPUTS.Hs = Hs;
    
%     save([NAME,'_inputs'],'NAME','INPUTS')
    
else % e.g. flag_loadCase == -1
    fprintf('\nRe-running the same case\n\n')
end

RESULTS = [];

for i = 1:n_mesh % n. times inc. is halved
    RESULTS(i).da_inc = da_tot*0.5^(i-1);
end


%% SOLVE CASE
for i = 1:n_mesh

% FIXED-LENGTH CRACK GROWTH SOLUTION USING A SINGLE TRIAL

da_h = zeros(n_tips,1); % fracture lengths
da_inc = RESULTS(i).da_inc; % increment length

n_grw = 0; % count number of crack increments until da_tot is achieved

while sum(da_h) < da_tot % number of crack increments
    
    Gs_h = Gs+Hs*da_h; % current energy release rates at crack tips
    da = zeros(n_tips,1)*da_inc; % trial solution: extend all tips uniformly
    q = abs((max(Gs_h)-Gs_h)/max(Gs_h)) < 0.01; % as in BudynZi2004
    
    if sum(q) > 1
        
        da(q) = CrackGrowthSol_MaxDeterminant_fixedLength(Gs_h(q),Hs(q,q),da_inc);
        
    else
        
        [~,q] = max(Gs_h); % time-explicit determination of critical crack tip
        da(q) = da_inc;
        
    end
    
    da_h = da_h + da;
    n_grw = n_grw + 1;
    
end

% APPROXIMATE IDEAL SOLUTION

da_tot = sum(da_h); % (corrected value)
da_inc = da_tot/n_grw; % (corrected value)
da_e = zeros(n_tips,1); % initialize exact

n_try = n_tips; % each trial corresponds to growth at each the crack tip

for j = 1:n_grw
    
    Gs_e = Gs+Hs*da_e; % current crack tip energy release rate
    
    ps_try = zeros(1,n_try);
    da_try = zeros(n_tips,n_try);
    
    for k = 1:n_try
        
        da = zeros(n_tips,1);
        da(k) = da_inc; % initial sol.
        
        [da_try(:,k),ps_try(k)] = CrackGrowthSol_GradientBasedOptim(Gs_e,Hs,da);
        
        % Note that any solution da in da_trials will satisfy:
        % sum(da(:,i)) == sum(da(:,j)) == da_inc for any i,j
        
    end
    
    [ps,k] = max(ps_try);
    da = da_try(:,k);
    da_e = da_e + da;
    
end

% COMPUTE DIFFERENCES

da_err = norm(da_h-da_e)/norm(da_e);
fprintf('Fracture solution discrepancy: %d\n',da_err)

if da_err < 2e-1
    break
end

Psi_e = da_e'*Gs+0.5*da_e'*Hs*da_e;
Psi_h = da_h'*Gs+0.5*da_h'*Hs*da_h;

RESULTS(i).da_e = da_e;
RESULTS(i).da_h = da_h;

RESULTS(i).da_err = da_err;
RESULTS(i).da_tot = da_tot;

RESULTS(i).Psi_e = Psi_e;
RESULTS(i).Psi_h = Psi_h;

end

if da_err > 2e-1
    break
end
end

fprintf('\nCASE SUMMARY:\n\n')
INPUTS
RESULTS
da_err
