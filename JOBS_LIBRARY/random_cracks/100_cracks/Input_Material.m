
%--------------------------------------------------------------------------
% Material data
%--------------------------------------------------------------------------

% 2D assumtions
problemType = 'PlaneStress'; % 'PlaneStress' or 'PlaneStrain'

% Units (plotting)
lengthUnits = 'mm';

% Material phases (to be consistent with mesh)
nPhase = 1;

% Material constants
E=[]; v=[];
E(1) = 1000; % MPa
v(1) = 0.3;

% Material toughness
K_crt = 1; % MPa.m^(1/2)

% unit conversion
switch lengthUnits
    case '\mum'
        K_crt = K_crt * 1e3; % MPa.m^0.5 == uN/um^2*um^0.5 * 1e3
    case 'mm'
        K_crt = K_crt * sqrt(1e3); % MPa.m^0.5 == N/mm^2*mm^0.5 * sqrt(1e3)
    case 'm'
        % okey.
    otherwise
        error('Length units?')
end

%--------------------------------------------------------------------------

