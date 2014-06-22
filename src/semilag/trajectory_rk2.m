function [xt,yt,zt] = trajectory_rk2(xx,yy,zz,fvel,ti,tf,n)
%TRAJECTORY_RK2 Computes the trajectory in a given velocity fields with
%               second order runge-kutta method.
if nargin < 7, n = 1; end;
global verbose;

dt = (tf - ti)/n;
xt = xx; yt = yy; zt = zz;
tt = ti;
for taustep=1:n
    [xt,yt,zt] = rk2(xt,yt,zt,fvel,tt,dt);
    if verbose
        plot(xt(:),yt(:),'gs','MarkerSize',8);  axis off; axis equal;
    end
    tt = tt + dt;
end
    %/* ************************************************** */
    function [xt,yt,zt] = rk2(xi,yi,zi,fvel,ti,dt)
        [u1,v1,w1] = fvel(ti,xi,yi,zi);
        k1x = dt*u1;
        k1y = dt*v1;
        k1z = dt*w1;
        tq = ti + 0.5*dt;
        xq = xi+0.5*k1x;
        yq = yi+0.5*k1y;
        zq = zi+0.5*k1z;
        [u2,v2,w2] = fvel(tq,xq,yq,zq);
        xt = xi + dt*u2;
        yt = yi + dt*v2;
        zt = zi + dt*w2;
    end
end

