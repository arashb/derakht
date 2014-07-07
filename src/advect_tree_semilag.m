%/* ************************************************** */
function [ctree_next] = advect_tree_semilag(ctree_next,fconc_interp,fvel_interp,t)
global resPerNode;

% PERFORM ONE STEP SEMI-LAGRANGIAN FOR EACH TREE LEAF
% TODO: remove tstep
tstep = 1;
cnext_leaves = ctree_next.leaves();
for lvcnt = 1:length(cnext_leaves)
    cnext_leaf  = cnext_leaves{lvcnt};    
    [xx,yy,zz,dx,dy,dz] = cnext_leaf.mesh(resPerNode);
    
    cnext_values = semilag_rk2(xx,yy,zz,fconc_interp,fvel_interp,t);
    
    cnext_leaf.data.dim           = 1;
    cnext_leaf.data.resolution    = resPerNode;
    cnext_leaf.data.values = cnext_values;
end
end