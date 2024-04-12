function [phi] = gradientflow(phi, maxiter, cname)
if nargin<2
    maxiter = 1e3;
end
global N c ep al epmckt2 dt dof kt kt2
tol = 1e-14;  res = Inf;  iter = 0;
%s = max(phi(:));
%fprintf('s = %.12e \n',s);
while res >= tol && iter < maxiter
    iter = iter + 1;
    phi1 = phi;
    phi2 = phi1 .^ 2 ;
    phi3 = phi2 .* phi1;
    ut = ifft2(reshape(phi1, N, N));  ut(1) = 0;

    nln = al * phi2 - phi3;
    nln = nln - mean(nln);
    %s1 = max(abs(nln(:)));
    %fprintf(' max nln %.3e',s1);
    nln = reshape(nln, N, N);
    if mod(iter, 5) == 0
        res = norm(ut + ifft2(nln) ./ epmckt2, 'fro') / N;
        ene = (c/2) * norm(kt.*ut, 'fro')^2 + ...
            (-(ep/2).*sum(phi2) - (al/3).*sum(phi3) + 0.25.*norm(phi2, 'fro')^2) / (N^2);
        drawcam(phi); 
        fprintf('%d: %e, E = %.12e\n', iter, res, ene );
    end
    
    ut = (ut + dt*ifft2(nln) )./( 1 - dt*epmckt2 );    % semi-implicit
    ut(1)=0;
%     mx = max(abs(ut(:)));
%     ut(abs(ut) < 0.0000001*mx)=0;
    phi1 = reshape(real(fft2(ut)), N^2, 1);
    phi = phi1(1:dof);
    
    %s = max(phi1(:));
    %fprintf('s = %.12e \n',s);
end
% fprintf('%d: %e, E = %.15e\n', iter, res, ene);
phi = phi1(1:dof);
% if nargin == 3
%     drawcam(phi, cname);
% else
%     drawcam(phi); 
%end

end