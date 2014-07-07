%/* ************************************************** */
function [ctree_next] = advect(c, ctree_next, ucells, vcells, t)
global DEBUG;
global resPerNode;
global INTERP_TYPE;

% MERGE VELOCITY TREES
num     = size(ucells,2);
um_tree = merge(ucells);
umcells = clone(um_tree,num);
vm_tree = merge(vcells);
vmcells = clone(vm_tree,num);

% INTERPOLATE VELOCITY VALUES ON THE MERGED TREES
for i =1:num
    tree_data.interp(ucells{i}, umcells{i});
    tree_data.interp(vcells{i}, vmcells{i});
end

um = tree_data.collapse(umcells);
vm = tree_data.collapse(vmcells);

if DEBUG
    n = length(ucells);
    figure('Name','MERGED VEL. TREES');
    for i=1:n
        subplot(3,n,i)
        ucells{i}.plottree;
        tree_data.plot_data(ucells{i});
        subplot(3,n,n+i)
        umcells{i}.plottree;
        tree_data.plot_data(umcells{i});
    end
    subplot(3,n,2*n+1);
    um_tree.plottree
    axis off; axis equal;
end

% PERFORM ONE STEP SEMI-LAGRANGIAN FOR EACH TREE LEAF
% TODO: remove tstep
tstep = 1;
cnext_leaves = ctree_next.leaves();
for lvcnt = 1:length(cnext_leaves)
    cnext_leaf  = cnext_leaves{lvcnt};    
    [xx,yy,zz,dx,dy,dz] = cnext_leaf.mesh(resPerNode);
    
%     fconc = @conc_exact;
%     fvel  = @vel_exact;
    fconc_interp = @conc_interp;
    fvel_interp  = @vel_interp;
    cnext_values = semilag_rk2(xx,yy,zz,fconc_interp,fvel_interp,t,tstep);
    
    cnext_leaf.data.dim           = 1;
    cnext_leaf.data.resolution    = resPerNode;
    cnext_leaf.data.values = cnext_values;
end
    
    %/* ************************************************** */
    function ci = conc_interp(tstep,xq,yq,zq)
        ci = tree_data.interp_points(c,xq,yq,zq);        
        ci = conc_out(ci,xq,yq,zq);
        
        function cq = conc_out(cq,xq,yq,zq)
            % OUTSIDE THE SIMULATION DOMAIN
            out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
            cq(out) = 0;
        end
    end

    %/* ************************************************** */
    function ci = conc_interp_deprecated(tstep,xt,yt,zt)
        ci = interp_conc_spatial(cnumsol,xx,yy,zz,tstep,xt,yt,zt,INTERP_TYPE,@conc_out);
        
        function cq = conc_out(cq,xq,yq,zq)                 
            % OUTSIDE THE CURRENT NODE
            [xmin,xmax,ymin,ymax] = c_leaf.corners();
            out = xq<xmin | xq>xmax  | yq<ymin | yq>ymax;% | zq<0 | zq>1;
            %ce = conc_exact(0,xq,yq,zq);
            ce = tree_data.interp_points(c,xq,yq,zq);
            cq(out) = ce(out);
            
            % OUTSIDE THE SIMULATION DOMAIN
            out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
            cq(out) = 0;
        end
    end

    %/* ************************************************** */
    function [uq,vq,wq] = vel_interp(tq,xq,yq,zq)
        uval = tree_data.interp_points(um,xq,yq,zq);
        vval = tree_data.interp_points(vm,xq,yq,zq);
        [uq,vq,wq,] = interp_vel_temporal(uval,vval,0,t,tq,INTERP_TYPE);
        [uq,vq,wq] = vel_out(uq,vq,wq,xq,yq,zq);
        
        function [uq,vq,wq] = vel_out(uq,vq,wq,xq,yq,zq)
            % OUTSIDE THE SIMULATION DOMAIN
            out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
            [ue, ve, we] = vel_rot(tq,xq,yq,zq,0.5,0.5,0.5);
            uq(out) = ue(out);
            vq(out) = ve(out);
            wq(out) = we(out);
        end
    end
    
    %/* ************************************************** */
    function [mt] = merge(treecells)
        mt = treecells{1};
        for counter=2:length(treecells)
            mt = qtree.merge(mt, treecells{counter});
        end
    end

    %/* ************************************************** */
    function [tree_clones] = clone(tree_src, num)
        tree_clones = cell(1,num);
        for counter=1:num
            tree_clones{counter} = qtree.clone(tree_src);
        end
    end
end