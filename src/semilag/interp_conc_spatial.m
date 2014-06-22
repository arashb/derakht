function [cq] = interp_conc_spatial(c,xx,yy,zz,tstep,xq,yq,zq,INTERP_TYPE,fout)
global dim;

cp = c(:,:,:,tstep);
switch dim
    case 3, cq = interp3(xx,yy,zz,cp,xq,yq,zq,INTERP_TYPE);
    case 2, cq = interp2(xx,yy,cp,xq,yq,INTERP_TYPE);
end
[cq] = fout(cq,xq,yq,zq);
end