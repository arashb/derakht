function [ cnew ] = semilag_rk2(xx,yy,zz,fconc,fvel,t,tstep)
%ADVECT_SL_RK2 Advect the c values one time step by using semi-lagrangian scheme
% 
VPREVTSTEP  = 1;
VCURTSTEP   = 2;
VNEXTSTEP   = 3;
V2NEXTSTEP  = 4;

ti = t(VNEXTSTEP);
tf = t(VCURTSTEP);
n  = 10;

[xt,yt,zt] = trajectory_rk2(xx,yy,zz,fvel,ti,tf,n);
cnew = fconc(tstep,xt,yt,zt);
end

