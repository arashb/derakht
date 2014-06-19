function [xt,yt,zt] = trajectory_rk2(xi,yi,zi,fvel,ti,tf,n)
%TRAJECTORY_RK2 Computes the trajectory in a given velocity fields with
%               second order runge-kutta method.
if nargin < 7, n = 1; end;
global verbose;

dt = (tf - ti)/n; 
xt = xi; yt = yi; zt = zi;
tt = ti;
for taustep=1:n 
    [xt,yt,zt] = rk2(xt,yt,zt,fvel,tt,dt);
    if verbose        
        plot(xt(:),yt(:),'gs','MarkerSize',8);  axis off; axis equal;
    end
    tt = tt + dt;
end
end

function [xt,yt,zt] = rk2(xi,yi,zi,fvel,ti,dt)
[u1,v1,w1] = fvel(ti,xi,yi,zi);
k1x = dt*u1;
k1y = dt*v1;
k1z = dt*w1;
[u2,v2,w2] = fvel( ti + 0.5*dt, xi+0.5*k1x, yi+0.5*k1y, zi+0.5*k1z);
xt = xi + dt*u2;
yt = yi + dt*v2;
zt = zi + dt*w2;
end