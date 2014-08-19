function [ csol ] = compute_analytical( xi, xf, ti, tf, xx, yy, zz, CF_TYPE )
%COMPUTE_ANALYTICAL  Compute/Load Analytical solution
global verbose;

fprintf('computing analytical solution:\n');
xc   = (xf + xi)/2;
yc   = xc;
zc   = xc;
csol = zeros(size(xx));
% compute the trajectory back in time from tf to ti
[xxorigin, yyorigin, zzorigin] = trajectory_ode45(xx, yy, zz, @vel, tf, ti);

if verbose
    figure
    MS='MarkerSize';
    plot(xx(:),yy(:),'o',MS,5); hold on;
    plot(xxorigin(:),yyorigin(:),'rx',MS,8);  axis off; axis equal;
    title('Traj. of Grid Points');
end

% init the concentration
switch CF_TYPE
    case 0
        csol = slotted_cylinder( xi, xf, xxorigin, yyorigin, zzorigin);
    case 1
        csol = gaussian(xxorigin, yyorigin, xc, yc, 0);
    case 2
        csol = cone(xxorigin, yyorigin);
end
assert(sum(sum(sum(isnan(csol)))) == 0,'NaN found in analytical solution data.');


    function [xt,yt,zt] = vel(t,x,y,z)
        [xt,yt,zt] = vel_rot(t,x,y,z,xc,yc,zc);
    end
end
