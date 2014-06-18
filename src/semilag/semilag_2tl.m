function [ cnew ] = semilag_2tl( interp_conc, xx, yy, zz, u, v, w, tstep, dt, vel_interp_method, conc_interp_method )

V2PREVTSTEP = 1;
VPREVTSTEP  = 2;
VCURTSTEP   = 3;
VNEXTSTEP   = 4;

% compute the velocity in midpoint (x - alpha/ 2, t + dt/2)
umidpoint = 1.5*u(:,:,:,VCURTSTEP) - 0.5*u(:,:,:,VPREVTSTEP);
vmidpoint = 1.5*v(:,:,:,VCURTSTEP) - 0.5*v(:,:,:,VPREVTSTEP);
wmidpoint = 1.5*w(:,:,:,VCURTSTEP) - 0.5*w(:,:,:,VPREVTSTEP);

% compute alpha: the distance the particle in semi-lagrangian is
% displaced in x in time dt

% inital value for fixed point iteration
% TODO: which velocity for initial value
xalpha = xx - u(:,:,:,VCURTSTEP)*dt;
yalpha = yy - v(:,:,:,VCURTSTEP)*dt;
zalpha = zz - w(:,:,:,VCURTSTEP)*dt;

fpicnt = 5;
for k=1:fpicnt
    xxtmp = xx - xalpha*0.5;
    yytmp = yy - yalpha*0.5;
    zztmp = zz - zalpha*0.5;
    
    ut = interp3(xx, yy, zz, umidpoint, xxtmp, yytmp, zztmp, vel_interp_method);
    vt = interp3(xx, yy, zz, vmidpoint, xxtmp, yytmp, zztmp, vel_interp_method);
    wt = interp3(xx, yy, zz, wmidpoint, xxtmp, yytmp, zztmp, vel_interp_method);
    
    out = xxtmp<0 | xxtmp>1  | yytmp<0 | yytmp>1 | zztmp<0 | zztmp>1;
    % TODO
    [ue, ve, we] = vel_rot(0,xxtmp,yytmp,zztmp,0.5,0.5,0.5);
    ut(out) = -ue(out);
    vt(out) = -ve(out);
    wt(out) = -we(out);
    
    xalpha = dt * ut;
    yalpha = dt * vt;
    zalpha = dt * wt;    

end


xt = xx - xalpha;
yt = yy - yalpha;
zt = zz - zalpha;

cnew = interp_conc(xt,yt,zt,tstep);
end

