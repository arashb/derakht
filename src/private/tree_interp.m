function tree_interp(src_tree, dst_tree)
% interpolate the values of grid points in the destination tree 
% from the values of correspoding nodes in the source tree
% TODO: this is a brute-force algorithm. optimize it.
% TODO: this algorithm works only if dst_tree is finer than source tree
global verbose;

dst_leaves = dst_tree.leaves();
src_leaves = src_tree.leaves();
resPerNode = src_leaves{1}.data.resolution;
for dst_lvcnt = 1:length(dst_leaves)
    dst_leaf = dst_leaves{dst_lvcnt};
    [dst_xxr,dst_yyr,dst_zzr,dst_dx,dst_dy,dst_dz] = dst_leaf.mesh(resPerNode);
    % TODO: find the correct dimension
    dst_leaf.data.dim = 1;
    dst_leaf.data.resolution = resPerNode;
    dst_leaf.data.values = zeros(size(dst_xxr));
    for src_lvcnt =1:length(src_leaves)
        src_leaf = src_leaves{src_lvcnt};
        indices = points_in_node(src_leaf, dst_xxr, dst_yyr);        
        if ~any(any(indices)), continue; end;
        xx = dst_xxr(indices);
        yy = dst_yyr(indices);
        if verbose
            fprintf('interpolate values of dst node %d from src node %d\n',dst_lvcnt, src_lvcnt);
        end
        [src_xxr,src_yyr,src_zzr,src_dx,src_dy,src_dz] = src_leaf.mesh(resPerNode);
        vv = interp2(src_xxr,src_yyr,src_leaf.data.values,xx,yy);
        dst_leaf.data.values(indices) = vv;    
    end
end

    %/* ************************************************** */
    function indices = points_in_node(node, xx, yy)
        % complication for points that lie right on
        % the boundaries
        [xmin,xmax,ymin,ymax]=corners(node);
        idx = find(xmin <= xx & xx < xmax);
        idy = find(ymin <= yy & yy < ymax);
        indices = intersect(idx, idy);
    end
end