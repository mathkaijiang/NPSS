% basic parameters 

% for LP model
global Q q0 q1
Q = 12; q0 = 1; q1 = 2 * cos(pi/Q);  

global c ep al dt
c = 1; ep = 0.05; al = 1.0;         dt = 1;   c1=-0.4;  c2=1.2;
%c = 24; ep =  1.2; al = 3;         dt = 0.1; c1=-0.4;  c2=1.4;
%c = 1; ep =  5e-6; al = 1/sqrt(2); dt = 1;   c1=-0.3;  c2=0.9;



global N L
%L = 30;  N = 512; 
%L = 82;  N = 1000;
L = 112;  N = 1024;


