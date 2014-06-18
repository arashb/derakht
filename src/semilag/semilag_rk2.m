function [ cnew ] = semilag_rk2(fconc_interp, fvel_interp, xx, yy, zz, t, tstep)
%ADVECT_SL_RK2 Advect the c values one time step by using semi-lagrangian scheme

V2PREVTSTEP = 1;
VPREVTSTEP  = 2;
VCURTSTEP   = 3;
VNEXTSTEP   = 4;

xi = xx;
yi = yy;
zi = zz;
ti = t(VCURTSTEP);
tf = t(VPREVTSTEP);
n = 10;

[xt,yt,zt] = trajectory_rk2(xi,yi,zi,fvel_interp,ti,tf,n);
cnew = fconc_interp(xt,yt,zt,tstep);

end

