function [ci] = interp_conc(c,xx,yy,zz,tstep,xt,yt,zt,INTERP_TYPE)
global dim;

cp = c(:,:,:,tstep);
switch dim
    case 3, ci = interp3(xx,yy,zz,cp,xt,yt,zt,INTERP_TYPE);
    case 2, ci = interp2(xx,yy,cp,xt,yt,INTERP_TYPE);
end
% TODO: fix this. gaussian function at time t is needed.
% ce = gaussian( xt, yt, zt, 0.5, 0.5, 0.5, 0);
out = xt<0 | xt>1  | yt<0 | yt>1 | zt<0 | zt>1;
ci(out) = 0;%ce(out);
end