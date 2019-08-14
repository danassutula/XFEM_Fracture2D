function vEnJmp = EnrJmp_Hvi(mGsCrd)

% Modified Heaviside function has a jump = 1 - (-1) = 2

jmp = 2;

vEnJmp = jmp(ones(size(mGsCrd,1),1),1);

end