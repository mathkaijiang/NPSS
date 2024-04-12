disp(['parameters: c=' num2str(c) '; ep=' num2str(ep) ';al=' num2str(al) ';'])
disp(['mesh: L=' num2str(L) ';N=' num2str(N) ';dt=' num2str(dt) ';'])

% set mesh  
global L N
global dof xgrids xs ys 
dof = N^2; 
xgrids = linspace(0, 2 * pi * L, N+1); xgrids(end)=[];
[ys, xs] = meshgrid([0:N/2-1 -N/2:-1]/L, [0:N/2-1 -N/2:-1]/L);
inpfp = @(p,q) (p'*q) / (N^2);      
nrmp = @(p) sqrt(inpfp(p,p));
rayleighq = @(x,v)inpfp(v,hv_camnew(x,v))/inpfp(v,v);
ldnm = @(NM)load(['./lp_camnew/data' num2str(c) '_' num2str(ep) '_' num2str(al) '/' num2str(L) '_' num2str(N)  '/S' num2str(L) '_' NM]);

% fourier coefficients
global kt ikt2 kt2
kt = ([0:N/2-1 -N/2:-1].^2 + [0:N/2-1 -N/2:-1]'.^2) ./ (L^2);
if q1 == 0
    kt = (q0^2 - kt); disp('LB model');
else
    kt = (q0^2 - kt) .* (q1^2 - kt); disp('LP model');
end

kt2 = kt.^2;

ikt2 = 1 ./ ( kt.^2 + abs(ep) +0.02); %ep/c

disp(['precond: p=' num2str(1/max(ikt2(:)))])
global epmckt2
epmckt2 = ep - c .* (kt.^2);

if ~exist(['./lp_camnew/data' num2str(c) '_' num2str(ep) '_' num2str(al) '/' num2str(L) '_' num2str(N) ], 'dir')
    mkdir(['./lp_camnew/data' num2str(c) '_' num2str(ep) '_' num2str(al) '/' num2str(L) '_' num2str(N) ]);
end
   