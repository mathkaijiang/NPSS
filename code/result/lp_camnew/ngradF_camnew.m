function F = ngradF_camnew(ut)
global epmckt2 al N dof
ut(1) = 0;
phi = real(fftn(ut));
phi = phi(:);
F = reshape( real(fft2(epmckt2.*ifft2(reshape(phi,N,N)))),N^2,1)...
    + (phi .^ 2) .* (al - phi);
F = F - mean(F);
F = F(1:dof);
end

