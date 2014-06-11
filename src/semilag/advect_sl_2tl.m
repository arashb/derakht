function [ cnew ] = advect_sl_2tl( c, xx, yy, zz, u, v, w, tstep, dt, vel_interp_method, conc_interp_method )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 8
    vel_interp_method =   'linear';
    conc_interp_method =   'linear';
else if  nargin == 9
        conc_interp_method =   'linear';
    end
end
VPREVTSTEP = 1;      % index for the vel. values of the previous time step
VCURTSTEP  = 2;      % index for the vel. values of the current time step
cnew = zeros(size(c(:,:,:,tstep)));
% compute the velocity in midpoint (x - alpha/ 2, t + dt/2)
umidpoint = 1.5*u(:,:,:,VCURTSTEP) - 0.5*u(:,:,:,VPREVTSTEP);
vmidpoint = 1.5*v(:,:,:,VCURTSTEP) - 0.5*v(:,:,:,VPREVTSTEP);
wmidpoint = 1.5*w(:,:,:,VCURTSTEP) - 0.5*w(:,:,:,VPREVTSTEP);
% compute alpha: the distance the particle in semi-lagrangian is
% displaced in x in time dt
% inital value for fixed point iteration
% TODO: which velocity for initial value
xalpha = xx(2:end-1,2:end-1,2:end-1) - u(2:end-1,2:end-1,2:end-1,VCURTSTEP)*dt;
yalpha = yy(2:end-1,2:end-1,2:end-1) - v(2:end-1,2:end-1,2:end-1,VCURTSTEP)*dt;
zalpha = zz(2:end-1,2:end-1,2:end-1) - w(2:end-1,2:end-1,2:end-1,VCURTSTEP)*dt;
xxtmp = xx(2:end-1,2:end-1,2:end-1) - xalpha;
yytmp = yy(2:end-1,2:end-1,2:end-1) - yalpha;
zztmp = zz(2:end-1,2:end-1,2:end-1) - zalpha;
fpicnt = 5;
for k=1:fpicnt
    xalpha = dt * interp3(xx, yy, zz, umidpoint, xxtmp, yytmp, zztmp, vel_interp_method);
    yalpha = dt * interp3(xx, yy, zz, vmidpoint, xxtmp, yytmp, zztmp, vel_interp_method);
    zalpha = dt * interp3(xx, yy, zz, wmidpoint, xxtmp, yytmp, zztmp, vel_interp_method);    
    xxtmp = xx(2:end-1,2:end-1,2:end-1) - xalpha;
    yytmp = yy(2:end-1,2:end-1,2:end-1) - yalpha;
    zztmp = zz(2:end-1,2:end-1,2:end-1) - zalpha;
end
% TODO: maybe you do not need the next 3 lines. CHECK IT!
xtmp = xx(2:end-1,2:end-1,2:end-1) - xalpha;
ytmp = yy(2:end-1,2:end-1,2:end-1) - yalpha;
ztmp = zz(2:end-1,2:end-1,2:end-1) - zalpha;
cnew(2:end-1,2:end-1,2:end-1) = interp3(xx,yy,zz,c(:,:,:,tstep),xtmp,ytmp,ztmp, conc_interp_method);
% pbc
cnew(1,:,:) = cnew(end-1,:,:);
cnew(end,:,:) = cnew(2,:,:);
cnew(:,1,:) = cnew(:,end-1,:);
cnew(:,end,:) = cnew(:,2,:);
cnew(:,:,1) = cnew(:,:,end-1);
cnew(:,:,end) = cnew(:,:,2);
end

