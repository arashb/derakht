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

%/* ************************************************** */
function [err] = refine_criterion_modified(qtree, func, resPerNode,t)
    global INTERP_TYPE;
    global VERBOSE;

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
    if VERBOSE, fprintf('-> REFINE_CRITERION: compute SL on grid centers!\n'); end;
    fce = func(t, xxcc, yycc, zzcc);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the interpolated values at the center of the local REGULAR grid cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the local grid points positions
    [xxg,yyg,zzg,dx,dy,dz] = qtree.mesh(resPerNode, INTERP_TYPE);

    if strcmp(INTERP_TYPE, 'CHEBYSHEV')
        if isempty(qtree.data),
            if VERBOSE, fprintf('-> REFINE_CRITERION: compute SL on grid points!\n'); end;
            fre = func(t,xxg,yyg,zzg);
            w = qdata.get_node_cheb_interpolant(qtree, fre, resPerNode);
        else
            w = qtree.data.values;
        end

        global CHEB_IMPL;
        if strcmp(CHEB_IMPL, 'IAS')
            [xmin,xmax,ymin,ymax] = qtree.corners;
            x = xxcc(1,:);
            y = yycc(:,1)';
            % scale the query points to -1 and 1
            xs = (x - xmin)*2/(xmax-xmin)-1.0;
            ys = (y - ymin)*2/(ymax-ymin)-1.0;
            fci = cheb.chebeval2(w,xs,ys);
        elseif strcmp(CHEB_IMPL, 'CHEBFUN')
            fci = w(xxcc, yycc);
        end
        %fci = w(xxcc, yycc);
    else
        if isempty(qtree.data),
            if VERBOSE, fprintf('-> REFINE_CRITERION: compute SL on grid points!\n'); end;
            fre = func(t,xxg,yyg,zzg);
        else
            fre = qtree.data.values;
        end
        fci = interp2(xxg, yyg, fre, xxcc, yycc,INTERP_TYPE);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute interpolation error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diff = fci - fce;
    err = max(abs(diff(:)));
end