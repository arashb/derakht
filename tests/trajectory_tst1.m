function trajectory_tst1()
clear; clear globals;
close all;
addpath('../src/common/');
addpath('../src/semilag/');
global verbose;
global gvfreq;

% RUN PARAMETERS
verbose         = false;
gvfreq          = 1;
ti              = 0.5;
tf              = 0;
xi              = 0;
xf              = 1;

% MAIN SCRIPT
fvel = @func1;

n = 2^5;
x = linspace(xi,xf,n+1);
[xr,yr,zr] = meshgrid(x,x,1);

% PLOT
MS='MarkerSize';
plot(xr(:),yr(:),'o',MS,5); hold on;
axis off; axis equal;
title('Traj. Points');

% exact trajectory
[xe,ye,ze] = trajectory_ode45(xr,yr,zr,fvel,ti,tf);
hold on;
plot(xe(:),ye(:),'rx',MS,8);  

% computed trajectory
tn = 10;
[xc,yc,zc] = trajectory_rk2(xr,yr,zr,fvel,ti,tf,tn);

hold on;
plot(xc(:),yc(:),'gs',MS,8);  axis off; axis equal;

end

function [xt,yt,zt] = func1(t,x,y,z)
    xc = 0.5;
    yc = xc;
    zc = 0;
    [xt,yt,zt] = vel_rot(t,x,y,z,xc,yc,zc);
end