function F = ngrad_camnew(phi)
global epmckt2 al N dof
F = reshape( real(fft2(epmckt2.*ifft2(reshape(phi,N,N)))),N^2,1)...
    + (phi .^ 2) .* (al - phi);
F = reshape(F,N,N);
Ff = ifftn(F);
Ff(1) =0;
F = real(fftn(Ff));
F = F(:);
%F = F - mean(F);
end