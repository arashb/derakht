function [uq, vq, wq] = interp_vel_spatial(xx,yy,zz,u,v,w,xq,yq,zq,INTERP_TYPE,fout)
%INTERP_VEL_SPATIAL Spatial interpolation of velocity values.
%
% ur, vr, wr are the velocity values at the grid points
% xq, yq, zq are the positions of queried velocities.
global dim;

switch dim
    case 3
        uq = interp3(xx,yy,zz,u,xq,yq,zq, INTERP_TYPE);
        vq = interp3(xx,yy,zz,v,xq,yq,zq, INTERP_TYPE);
        wq = interp3(xx,yy,zz,w,xq,yq,zq, INTERP_TYPE);
    case 2
        uq = interp2(xx,yy,u,xq,yq, INTERP_TYPE);
        vq = interp2(xx,yy,v,xq,yq, INTERP_TYPE);
        wq = interp2(xx,yy,w,xq,yq, INTERP_TYPE);
end
[uq,vq,wq] = fout(uq,vq,wq,xq,yq,zq);
end