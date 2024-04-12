clc
clear
addpath lp_camnew
close all
set_para;

initialize_cam;
%choice = 8; inidegree =0 ;
%[x, cname] = guesses(choice, inidegree);


%[x] = gradientflow(x, 1000, cname);
%V = load('V1_GOSD.txt');
%V = reshape(V,dof,5);
%drawcam(V(:,1));


x = load('x1.txt');
drawcam(x);
drawnow
% 
%initeigs;
%giveop;