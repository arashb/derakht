function [cq] = interp_conc_spatial(c,xx,yy,zz,xq,yq,zq,INTERP_TYPE,fout)
global dim;

switch dim
    case 3, cq = interp3(xx,yy,zz,c,xq,yq,zq,INTERP_TYPE);
    case 2, cq = interp2(xx,yy,c,xq,yq,INTERP_TYPE);
end
[cq] = fout(cq,xq,yq,zq);
end