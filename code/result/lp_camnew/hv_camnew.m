function hv = hv_camnew(phi,v)
global epmckt2 al N dof
hv = -reshape(real(fftn(epmckt2.*ifftn(reshape(v,N,N)))), N^2,1)...
    + (3*phi - 2*al) .* phi .* v;

%fprintf(" mean: %.3e",mean(hv(:)));
hv = hv - mean(hv(:));

end

