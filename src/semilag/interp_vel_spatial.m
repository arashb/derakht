function [uq, vq, wq] = interp_vel_spatial(xx,yy,zz,u,v,w,xq,yq,zq,INTERP_TYPE)
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

out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
% TODO
[ue, ve, we] = vel_rot(0,xq,yq,zq,0.5,0.5,0.5);
uq(out) = -ue(out);
vq(out) = -ve(out);
wq(out) = -we(out);
end