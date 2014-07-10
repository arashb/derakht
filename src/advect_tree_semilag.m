%/* ************************************************** */
function cnext = advect_tree_semilag(c,fvelx,fvely,t,fdo_refine,fconc_exact,fvel_exact)
global resPerNode;
global verbose;
global INTERP_TYPE;

VPREVTSTEP  = 1;
VCURTSTEP   = 2;
VNEXTSTEP   = 3;
V2NEXTSTEP  = 4;

% VELOCITY (TIME-DEPENDENT)
nt = length(t);
ucells = cell(1,nt);
vcells = cell(1,nt);

for tcnt = 1:nt
    utmptree = qtree;
    utmptree.insert_function(fvelx,fdo_refine,t(tcnt));
    tree_data.init_data(utmptree,fvelx,resPerNode,t(tcnt));
    ucells{tcnt} = utmptree;
    
    vtmptree = qtree;
    vtmptree.insert_function(fvely,fdo_refine,t(tcnt));
    tree_data.init_data(vtmptree,fvely,resPerNode,t(tcnt));
    vcells{tcnt} = vtmptree;
end

% MERGE VELOCITY TREES
fprintf('-> merge velocity trees\n');
num     = size(ucells,2);
um_tree = merge(ucells);
umcells = clone(um_tree,num);
vm_tree = merge(vcells);
vmcells = clone(vm_tree,num);

% INTERPOLATE VELOCITY VALUES ON THE MERGED TREES
fprintf('-> interpolate velocity values on the merged tree\n');
for i =1:num
    tree_data.interp(ucells{i}, umcells{i});
    tree_data.interp(vcells{i}, vmcells{i});
end

um = tree_data.collapse(umcells);
vm = tree_data.collapse(vmcells);

% ADVECT
fprintf('-> performing one step semi-lagrangian\n');
%fconc_interp    = @conc_exact;
fconc_interp    = @conc_interp;
%fvel_interp     = @vel_exact;
fvel_interp     = @vel_interp;
%fsemilag        = fconc_exact;
fsemilag        = @semilag;

tic
% FIRST METHOD: CONSTRUCT THE NEXT TIME STEP TREE FROM SCRATCH WITH
%               SEMILAG SOLVER AS REFINEMENT FUNCTION
% cnext = qtree;
% cnext.insert_function(fsemilag,fdo_refine);
% tree_data.init_data(cnext,fsemilag,resPerNode);

% SECOND METHOD: USE THE PREVIOUS TIME STEP TREE AS STARTING POINT
%                REFINE/COARSEN WHENEVER IS NEEDED
cnext = qtree.clone(c);
update_tree(cnext,fdo_refine);
toc

if verbose
    plotvel();
end

    function val = update_tree(node, fvisit)
        val = true;
        kidsval = true;
        if ~node.isleaf
            for k=1:4
                kidval = update_tree(node.kids{k}, fvisit);
                kidsval = kidsval & kidval;
            end
        end
        [refine_node, values] = fvisit(node,fsemilag,t(VNEXTSTEP));
        % REFINE THE NODE
        if refine_node & isempty(node.kids())
            if verbose,
                mid = morton_id;
                id = mid.id(node.level,node.anchor);
                fprintf('--> refine node: ')
                mid.print(id)
                fprintf(' level %2d: anchor:[%1.4f %1.4f] \n', ...
                    node.level,node.anchor(1),node.anchor(2));
            end
            refine(node);
            for kcnt=1:4, update_tree(node.kids{kcnt},fvisit); end;
            val = false;
            return
        end
        % COARSEN THE NODE
        if ~refine_node & kidsval & ~isempty(node.kids())
            if verbose,
                mid = morton_id;
                id = mid.id(node.level,node.anchor);
                fprintf('--> coarsen node: ')
                mid.print(id)
                fprintf(' level %2d: anchor:[%1.4f %1.4f] \n', ...
                    node.level,node.anchor(1),node.anchor(2));
            end
            coarsen(node)
            set_node_values(node, values);
            return;
        end
        % KEEP THE NODE AS IT IS
        set_node_values(node, values);
        
        function set_node_values(node, values)
            if verbose,
                mid = morton_id;
                id = mid.id(node.level,node.anchor);
                fprintf('--> set semilag values for node: ')
                mid.print(id)
                fprintf(' level %2d: anchor:[%1.4f %1.4f] \n', ...
                    node.level,node.anchor(1),node.anchor(2));
            end
            node.data.dim           = 1;
            node.data.resolution    = resPerNode;
            node.data.values = values;
        end
    end

    function refine(node)
        node.data = [];
        node.create_kids();
    end

    function coarsen(node)
        node.kids = [];
        node.isleaf = true;
    end

    function val = semilag(tdummy,x,y,z)
        val = semilag_rk2(x,y,z,fconc_interp,fvel_interp,t);
    end

%/* ************************************************** */
    function ci = conc_interp(tq,xq,yq,zq)
        ci = tree_data.interp_points(c,xq,yq,zq);
        ci = conc_out(ci,xq,yq,zq);
        
        function cq = conc_out(cq,xq,yq,zq)
            % OUTSIDE THE SIMULATION DOMAIN
            ce = fconc_exact(tq,xq,yq,zq);
            out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
            cq(out) = ce(out);
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
            [ue, ve, we] = fvel_exact(tq,xq,vq,zq);
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

    function plotvel()
        figure('Name','SEMI-LAG QUAD-TREES');        
        subplot(3,4,2)
        c.plottree;
        tree_data.plot_data(c);
        title('c(t)');
        
        subplot(3,4,3)
        cnext.plottree;
        tree_data.plot_data(cnext);
        title('c(t+dt)');
        
        subplot(3,4,5)
        ucells{1}.plottree;
        tree_data.plot_data(ucells{1});
        title('u(t(n-1))');
        
        subplot(3,4,6)
        ucells{2}.plottree;
        tree_data.plot_data(ucells{2});
        title('u(t(n))');
        
        subplot(3,4,7)
        ucells{3}.plottree;
        tree_data.plot_data(ucells{3});
        title('u(t(n+1))');
        
        subplot(3,4,8)
        ucells{4}.plottree;
        tree_data.plot_data(ucells{4});
        title('u(t(n+2))');
        
        subplot(3,4,9)
        vcells{1}.plottree;
        tree_data.plot_data(vcells{1});
        title('v(t(n-1))');
        
        subplot(3,4,10)
        vcells{2}.plottree;
        tree_data.plot_data(vcells{2});
        title('v(t(n))');
        
        subplot(3,4,11)
        vcells{3}.plottree;
        tree_data.plot_data(vcells{3});
        title('v(t(n+1))');
        
        subplot(3,4,12)
        vcells{4}.plottree;
        tree_data.plot_data(vcells{4});
        title('v(t(n+2))');
    end
end