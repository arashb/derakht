%/* ************************************************** */
function advect_tree_semilag(cnext,fconc_interp,fvel_interp,t)
global resPerNode;
global verbose;

% PERFORM ONE STEP SEMI-LAGRANGIAN FOR EACH TREE LEAF
cnext_leaves = cnext.leaves();
for lvcnt = 1:length(cnext_leaves)
    cnext_leaf  = cnext_leaves{lvcnt};
    if verbose,
        mid = morton_id;
        id = mid.id(cnext_leaf.level,cnext_leaf.anchor);
        fprintf('--> compute semilag for leaf: ')
        mid.print(id)
        fprintf(' level %2d: anchor:[%1.4f %1.4f] \n', ...
            cnext_leaf.level,cnext_leaf.anchor(1),cnext_leaf.anchor(2));
    end
    [xx,yy,zz,dx,dy,dz] = cnext_leaf.mesh(resPerNode);
    
    cnext_values = semilag_rk2(xx,yy,zz,fconc_interp,fvel_interp,t);
    
    cnext_leaf.data.dim           = 1;
    cnext_leaf.data.resolution    = resPerNode;
    cnext_leaf.data.values = cnext_values;
end
end