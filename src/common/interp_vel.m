function [uq,vq,wq] = interp_vel(xx,yy,zz,u,v,w,t,tq,xq,yq,zq,INTERP_TYPE)
[ut, vt, wt] = interp_vel_temporal(u,v,w,t,tq,INTERP_TYPE);
[uq, vq, wq] = interp_vel_spatial(xx,yy,zz,ut,vt,wt,xq,yq,zq,INTERP_TYPE);
end