%% Operador (I-R_w)
function result = GV_LS(f, Kj_11,Kj_12,Kj_21,Kj_22,V_11,V_12,V_21,V_22, M)
% Reshape w to vertical vector
f = f(:);
f1 = f(1:2^(2*M));
f2 = f(2^(2*M)+1:2^(2*M+1));

% Temporary variable for w given on the square grid
f1 = reshape(f1, 2^M, 2^M);
f2 = reshape(f2, 2^M, 2^M);

ff1 = fftshift(fftn(fftshift(f1.*V_11+f2.*V_12)));
ff2 = fftshift(fftn(fftshift(f1.*V_21+f2.*V_22)));

% Multiply by the Fourier transform of the Green function and invert the FT
result1 = fftshift(ifftn(fftshift(Kj_11.*ff1 + Kj_12.*ff2)));
result2 = fftshift(ifftn(fftshift(Kj_21.*ff1 + Kj_22.*ff2)));
result=[result1(:);result2(:)];

% Final result needs the identity operator
result = f - result;
end
