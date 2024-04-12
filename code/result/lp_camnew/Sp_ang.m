function a = Sp_ang(V0,v)
global N 
inpfp = @(p,q) ([p;-sum(p)]'*[q;-sum(q)]) / (N^2);      
nrmp = @(p) sqrt(inpfp(p,p));
v1 = v - V0*inpfp(V0,v);
v1 = nrmp(v1);
a = asin(v1)*180/pi;
end