%/* ************************************************** */
function [val, funcval] = tree_do_refine(qtree, func, maxErrPerNode, maxLevel, resPerNode, t)
[err, funcval] = refine_criterion(qtree, func, resPerNode,t);

if qtree.level == maxLevel || err <=  maxErrPerNode
    val = false;
    return;
end

val = true;
end

function [err, fre]= refine_criterion(qtree, func, resPerNode,t)
global INTERP_TYPE;

[xxr,yyr,zzr,dx,dy,dz] = qtree.mesh(resPerNode);

% compute the function values on the local grid points
fre = func(t,xxr,yyr,zzr);

% compute the center of the local grid cells
xxc = xxr+dx/2;
yyc = yyr+dy/2;
zzc = zzr+dz/2;
xxcc = xxc(1:end-1,1:end-1);
yycc = yyc(1:end-1,1:end-1);
zzcc = zzc(1:end-1,1:end-1);

% compute the exact values on the centers
fce = func(t, xxcc, yycc, zzcc);

% interpolate the function values on the center points
fci = interp2(xxr, yyr, fre, xxcc, yycc,INTERP_TYPE);

% compute interpolation error
diff = fci - fce;
err = max(max(abs(diff)));
end