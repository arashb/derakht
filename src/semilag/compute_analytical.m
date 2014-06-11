function [ csol ] = compute_analytical( xi, xf, ti, tf, xx, yy, zz, CF_TYPE )
%COMPUTE_ANALYTICAL  Compute/Load Analytical solution
%
% In case that analytical solution does not exist
% computes it, otherwise the existing data will be loaded.
fprintf('computing analytical solution:\n');
xc = (xf + xi) / 2;
yc = xc;
zc = xc;
csol = zeros(size(xx));
% compute the trajectory back in time from tf to ti
[xxorigin, yyorigin, zzorigin] = trajectory(xx, yy, zz, @vel_rot, ti, tf);

% figure
% MS='MarkerSize';
% plot(xx(:),yy(:),'o',MS,5); hold on;
% plot(xxorigin(:),yyorigin(:),'rx',MS,8);  axis off; axis equal;
% title('Traj. Points of ODE45');

% init the concentration
switch CF_TYPE
    case 0
        csol = slotted_cylinder( xi, xf, xxorigin, yyorigin, zzorigin, dx);
    case 1
        csol = gaussian(xxorigin, yyorigin, xc, yc, 0);
end
assert(sum(sum(sum(isnan(csol)))) == 0,'NaN found in analytical solution data.');
end

