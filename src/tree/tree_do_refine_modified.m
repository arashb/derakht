%/* ************************************************** */
function [val] = tree_do_refine_modified(qnode, func, maxErrPerNode, maxLevel, resPerNode, t)
if qnode.level >= maxLevel
    val = false;
    return;
end

err = refine_criterion_modified(qnode, func, resPerNode,t);

if err <=  maxErrPerNode
    val = false;
    return;
end

val = true;
end

function [err]= refine_criterion_modified(qtree, func, resPerNode,t)
global INTERP_TYPE;

[xxr,yyr,zzr,dx,dy,dz] = qtree.mesh(resPerNode);

% get local grid points values
fre = qtree.data.values;

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