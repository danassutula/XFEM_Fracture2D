function [D,alpha,kappa] = MtxConstit(n,E,v,probtype)

D       =  cell(n,1);
alpha   = zeros(n,1);
kappa   = zeros(n,1);

switch probtype
    case 'PlaneStress'
        
        for i = 1:n
            
            % stress - strain relationship
            D{i} = [    1,  v(i),          0;   ...
                     v(i),     1,          0;   ...
                        0,     0, (1-v(i))/2    ].*(E(i)/(1-v(i)^2));
            
            % fracture parameters
            kappa(i) = (3-v(i))/(1+v(i));
            alpha(i) = 1;
            
        end
        
    case 'PlaneStrain'
        
        for i = 1:n
            
            % stress - strain relationship
            D{i} = [ 1-v(i),    v(i),       0;  ...
                       v(i),  1-v(i),       0;  ...
                          0,       0, 0.5-v(i) ].*(E(i)/(1+v(i))/(1-2*v(i)));
            
            % fracture parameters
            kappa(i) = 3-4*v(i);
            alpha(i) = 1-v(i)^2;
            
        end
        
    otherwise
        error('Problem type: ''PlaneStress'' or ''PlaneStrain'' ?')
        
end
end
