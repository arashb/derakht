%/* ************************************************** */
function [val, funcval] = tree_do_refine(qtree, func, maxErrPerNode, maxLevel, resPerNode, t)
    [err, funcval] = refine_criterion(qtree, func, resPerNode,t);

    if qtree.level == maxLevel || err <=  maxErrPerNode
        val = false;
        return;
    end

    val = true;
end

function [err, fre] = refine_criterion(qtree, func, resPerNode,t)
    global INTERP_TYPE;
    global verbose;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the exact values at the center of the local REGULAR grid cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the center of the local REGULAR grid cells
    [xxr,yyr,zzr,dx,dy,dz] = qtree.mesh(resPerNode);
    xxc = xxr+dx/2;
    yyc = yyr+dy/2;
    zzc = zzr+dz/2;
    xxcc = xxc(1:end-1,1:end-1);
    yycc = yyc(1:end-1,1:end-1);
    zzcc = zzc(1:end-1,1:end-1);

    % compute the exact values on the centers
    if verbose, fprintf('-> REFINE_CRITERION: compute SL on grid centers!\n'); end;
    fce = func(t, xxcc, yycc, zzcc);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the interpolated values at the center of the local REGULAR grid cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the local grid points positions
    [xxg,yyg,zzg,dx,dy,dz] = qtree.mesh(resPerNode, INTERP_TYPE);

    % compute the function values on the local grid points
    if verbose, fprintf('-> REFINE_CRITERION: compute SL on grid points!\n'); end;
    fre = func(t,xxg,yyg,zzg);

    % interpolate the function values on the center points
    if strcmp(INTERP_TYPE, 'CHEBYSHEV')
        [xmin xmax ymin ymax] = qtree.corners;
        w = chebfun2(fre, [xmin xmax ymin ...
                            ymax]);
        fci = w(xxcc, yycc);
    else
        fci = interp2(xxg, yyg, fre, xxcc, yycc,INTERP_TYPE);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute interpolation error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diff = fci - fce;
    err = max(abs(diff(:)));
end