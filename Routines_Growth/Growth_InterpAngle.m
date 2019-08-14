function b = Growth_InterpAngle(KK,KK_ref,b_ref)

% KK = K2/K1;
s  = sign(KK);
KK = s*KK;

% since monotonically decreasing
i = find(KK_ref > KK,1,'last');

if i > 1
    
    % interpolate
    b = (b_ref(i)*(KK_ref(i+1)-KK) + b_ref(i+1)*...
        (KK-KK_ref(i)))/(KK_ref(i+1)-KK_ref(i));
    
else % KK_ref(1) == inf
    
    % max angle
    b = b_ref(1);
    
end

b = b*s;

end