function [e0_xpt,s0_xpt] = Input_BC_preload(x,input_args)

%--------------------------------------------------------------------------
% Pre-stress/strain at discrete layers
%--------------------------------------------------------------------------

% Physically meaningful to specify only one: either stress or strain;
% therefore, return the other as empty or zeros. The return size is: [],
% [3,1] or [3,n], where n = size(x,1).

% x = element-center coordinates where pre-load is computed, size(x)=[n,2]

% EXAMPLE
% job_preload_args{1} = job_BC_preload; % pre-strain tensor
% job_preloadr_args{2} = 0.5*(job_domain_yi(4)+job_domain_yi(5));
% job_preload_args{3} = job_domain_yi(5)-job_domain_yi(4)];

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Simple case: a single pre-stress/strain for the whole domain
%--------------------------------------------------------------------------

if length(input_args) == 1
    e0_xpt = repmat(input_args{1}(:),1,size(x,1));
    s0_xpt = [];
    return
end

%--------------------------------------------------------------------------
% Sort input arguments (user-adapted)
%--------------------------------------------------------------------------

% layer pre-strain: [exx_i,eyy_i,exy_i; ...]
e0_lay = input_args{1};

% layer pre-stress: [sxx_i,syy_i,sxy_i; ...]
s0_lay = [];

% position of pre-load layers
depth = input_args{2};

% thickness of pre-load layers
thick = input_args{3};

%--------------------------------------------------------------------------

e0_xpt = [];
s0_xpt = [];

n = size(x,1);
N = length(depth);

if length(thick)~=N || (size(e0_lay,1)~=N && size(s0_lay,1)~=N)
    error('Inconsistent size between layer depth, thickness and pre-load.');
end

if ~isempty(e0_lay)
    e0_xpt = zeros(3,n);
    for i = 1:N
        
        q = find(x(:,2) > depth(i)-0.5*thick(i) ...
            & x(:,2) < depth(i)+0.5*thick(i));
        
        e0_xpt(1,q) = e0_lay(i,1);
        e0_xpt(2,q) = e0_lay(i,2);
        e0_xpt(3,q) = e0_lay(i,3);
        
    end
else % ~isempty(s0_lay)
    s0_xpt = zeros(3,n);
    for i = 1:N
        
        q = find(x(:,2) > depth(i)-0.5*thick(i) ...
            & x(:,2) < depth(i)+0.5*thick(i));
        
        s0_xpt(1,q) = s0_lay(i,1);
        s0_xpt(2,q) = s0_lay(i,2);
        s0_xpt(3,q) = s0_lay(i,3);
        
    end
end
end
