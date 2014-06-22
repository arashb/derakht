function [ cnew ] = semilag_rk2(xx,yy,zz,fconc,fvel,t,tstep)
%ADVECT_SL_RK2 Advect the c values one time step by using semi-lagrangian scheme
% 
V2PREVTSTEP = 1;
VPREVTSTEP  = 2;
VCURTSTEP   = 3;
VNEXTSTEP   = 4;

ti = t(VCURTSTEP);
tf = t(VPREVTSTEP);
n  = 10;

[xt,yt,zt] = trajectory_rk2(xx,yy,zz,fvel,ti,tf,n);
cnew = fconc(tstep,xt,yt,zt);
end

