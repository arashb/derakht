function [ cnew ] = semilag_rk2(c,xx,yy,zz,fconc_interp,u,v,w,t,fvel_interp,tstep,INTERP_TYPE)
%ADVECT_SL_RK2 Advect the c values one time step by using semi-lagrangian scheme

V2PREVTSTEP = 1;
VPREVTSTEP  = 2;
VCURTSTEP   = 3;
VNEXTSTEP   = 4;

ti = t(VCURTSTEP);
tf = t(VPREVTSTEP);
n  = 10;

[xt,yt,zt] = trajectory_rk2(xx,yy,zz,u,v,w,t,fvel_interp,ti,tf,n,INTERP_TYPE);
cnew = fconc_interp(c,xx,yy,zz,tstep,xt,yt,zt,INTERP_TYPE);
end

