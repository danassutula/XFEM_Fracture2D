function [da,ps] = CrackGrowthSol_MaxDeterminant_fixedLength(Gs,Hs,da_inc)
%
% a = assumed uniform extension of crack tips
% D = subdeterminant of Hs corresponding to a
% Hs = rates of the energy release rates

N = size(Hs,1); % n. competing crack tips
D = -inf; % initialise subdeterminant of Hs
a = false(1,N); % initialise crack extensions

for i = 1:2^N-1 % test all cases
    
    % get logical indices into Hs
    ai = logical(sprintf(['%0',num2str(N),'s'],dec2bin(i)) - '0');
    
    % get subdeterminant
    Di = det(Hs(ai,ai));
    
    if D < Di
        D = Di;
        a = ai;
    end
    
end

da = a(:)*da_inc;
ps = da'*(Gs+Hs*da*0.5);
