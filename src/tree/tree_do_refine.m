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
    % get the local grid points positions
    [xxg,yyg,zzg,dx,dy,dz] = qtree.mesh(resPerNode, INTERP_TYPE);

    % compute the function values on the local grid points
    if VERBOSE, fprintf('-> REFINE_CRITERION: compute SL on grid points!\n'); end;
    fre = func(t,xxg,yyg,zzg);

    % interpolate the function values on the center points
    if strcmp(INTERP_TYPE, 'CHEBYSHEV')
        global CHEB_IMPL;
        w = qdata.get_node_cheb_interpolant(qtree, fre, resPerNode);
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
    else
        fci = interp2(xxg, yyg, fre, xxcc, yycc,INTERP_TYPE);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute interpolation error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diff = fci - fce;
    [err, indx] = max(abs(diff(:)));
    [i_row, i_col] = ind2sub(size(diff),indx);
    % global MAX_ERROR_PER_NODE;
    % if err > maxErrorPerNode
    %     fci
    %     fce
    %     fprintf('COORD: (%2d, %2d)\n',         i_row, i_col);
    %     fprintf('INTRP: %1.4f\n',         fci(i_row, i_col));
    %     fprintf('EXACT: %1.4f\n',         fce(i_row, i_col));
    %     fprintf('ERROR: %1.4f\n',         err);
    %     fprintf('NODE: level %2d: anchor:[%1.4f %1.4f] \n', ...
    %             qtree.level,qtree.anchor(1),qtree.anchor(2));
    % end
end