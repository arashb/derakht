function main()
clear; clear globals; % constants  and preamble
addpath semilag
addpath tree
addpath common

global verbose;
global gvfreq;
global dim;
global DEBUG;
global resPerNode;
global maxErrorPerNode;
global maxLevel;
global INTERP_TYPE;

% RUN PARAMETERS
maxErrorPerNode = 0.001;        % Error per box
maxLevel        = 20;           % Maximum tree depth
resPerNode      = 15;           % Resolution per Node
verbose         = true;
gvfreq          = 1;
dim             = 2;
DEBUG           = false;
INTERP_TYPE     = 'cubic';

% MAIN SCRIPT
fconc_exact   = @conc_exact_2;
fvelx_exact   = @velx_exact;
fvely_exact   = @vely_exact;

% CONSTRUCT AND INIT THE TREES VALUES
fprintf('-> init the trees\n');

% CONCENTRATION
c = qtree;
tinit = 0;
c.insert_function(fconc_exact,@do_refine,tinit);
tree_data.init_data(c,fconc_exact,resPerNode,tinit);

% SIMULATION PARAMETERS BASED ON INITIAL CCONC. TREE
cdepth      = c.find_depth();
cwidth      = 1/2^cdepth;
dx          = cwidth/resPerNode;
cfl         = 100;
om          = 1;
dt          = cfl*dx/om;
VPREVTSTEP  = 1;
VCURTSTEP   = 2;
VNEXTSTEP   = 3;
V2NEXTSTEP  = 4;
t           = [-dt 0 dt 2*dt];

fprintf('--> tree depth: %d\n',cdepth);
fprintf('--> cfl: %d\n',cfl);
fprintf('--> dt: %d\n',dt);

% VELOCITY (TIME-DEPENDENT)
nt = length(t);
ucells = cell(1,nt);
vcells = cell(1,nt);

for tcnt = 1:nt
    utmptree = qtree;
    utmptree.insert_function(fvelx_exact,@do_refine,t(tcnt));
    tree_data.init_data(utmptree,fvelx_exact,resPerNode,t(tcnt));
    ucells{tcnt} = utmptree;
    
    vtmptree = qtree;
    vtmptree.insert_function(fvely_exact,@do_refine,t(tcnt));
    tree_data.init_data(vtmptree,fvely_exact,resPerNode,t(tcnt));
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
%fvel_interp     = @vel_exact;
fconc_interp    = @conc_interp;
fvel_interp     = @vel_interp;
%fsemilag        = fconc_exact;
fsemilag        = @semilag;

tic

% FIRST METHOD: CONSTRUCT THE NEXT TIME STEP TREE FROM SCRATCH WITH 
%               SEMILAG SOLVER AS REFINEMENT FUNCTION
% cnext = qtree;
% cnext.insert_function(fsemilag,@do_refine);
% tree_data.init_data(cnext,fsemilag,resPerNode);

% SECOND METHOD: USE THE PREVIOUS TIME STEP TREE AS STARTING POINT
%                REFINE/COARSEN WHENEVER IS NEEDED
cnext = qtree.clone(c);
update_tree(cnext,@do_refine);
toc

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

% PLOT THE RESULTS
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

f = figure('Name','SEMI-LAG ADVECTION');
subplot(1,2,1);
c.plottree(0.5);
tree_data.plot_data(c)
colorbar;
title('c(t)');

subplot(1,2,2);
cnext.plottree(0.5);
tree_data.plot_data(cnext);
colorbar;
title('c(t+dt)');

s_fig_name = ['results_', datestr(now)];
saveas(f,s_fig_name,'pdf');
    

    function val = semilag(tdummy,x,y,z)
        val = semilag_rk2(x,y,z,fconc_interp,fvel_interp,t);
    end

    %/* ************************************************** */
    function ci = conc_interp(t,xq,yq,zq)
        ci = tree_data.interp_points(c,xq,yq,zq);        
        ci = conc_out(ci,xq,yq,zq);
        
        function cq = conc_out(cq,xq,yq,zq)
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
    function [val, funcval] = do_refine(qtree,func,t)
        [val, funcval] = tree_do_refine(qtree,func,maxErrorPerNode,maxLevel,resPerNode,t);
    end

    %/* ************************************************** */
    function val = do_coarsen(qtree,func,t)
        tmpval = tree_do_refine(qtree,func,maxErrorPerNode,maxLevel,resPerNode,t);
        val = ~tmpval;
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

%/* ************************************************** */
function value = conc_exact_1(t,x,y,z)
xc = 0.5;
yc = 0.5;
om = 1;
theta = -om*t;
sigmax = 0.05;
sigmay = 0.12;
value = gaussian(x,y,xc,yc,theta, sigmax, sigmay);
end

%/* ************************************************** */
function value = conc_exact_2(t,x,y,z)
xc = 0.5;
yc = 0.5;
xci = 0.75;
yci = 0.25;
om = 1;

[alpha,RHO] = cart2pol(xci-xc,yci-xc);
alphat = alpha + t*om;

[xct,yct] = pol2cart(alphat,RHO);
xct = xct + xc;
yct = yct + yc;

theta = 0;
sigmax = 0.05;
sigmay = 0.05;
value = gaussian(x,y,xct,yct,theta,sigmax,sigmay);
end

%/* ************************************************** */
function value = conc_exact_22(t,x,y,z)
xc = 0.5;
yc = 0.5;
theta = 0;
sigmax = 0.05;
sigmay = 0.05;
om = 1;

xc1 = 0.75;
yc1 = 0.75;
[alpha1,RHO1] = cart2pol(xc1-xc,yc1-xc);
alphat1 = alpha1 + t*om;
[xct,yct] = pol2cart(alphat1,RHO1);
xct = xct + xc;
yct = yct + yc;
value1 = gaussian(x,y,xct,yct,theta,sigmax,sigmay);

xc2 = 0.25;
yc2 = 0.25;
[alpha2,RHO2] = cart2pol(xc2-xc,yc2-xc);
alphat2 = alpha2 + t*om;
[xct,yct] = pol2cart(alphat2,RHO2);
xct = xct + xc;
yct = yct + yc;
value2 = gaussian(x,y,xct,yct,theta,sigmax,sigmay);

value = value1 + value2;

end

%/* ************************************************** */
function value = conc_exact_3(t,x,y,z)
xi = 0;
xf = 1;
value = slotted_cylinder( xi, xf, x, y, z);
end

%/* ************************************************** */
function [u,v,w] = vel_exact(t,x,y,z)
    xc = 0.5;
    yc = 0.5;
    zc = 0.5;
    [u,v,w] = vel_rot(t,x,y,z,xc,yc,zc);
end

%/* ************************************************** */
function u = velx_exact(t,x,y,z)
    [u,v,w] = vel_exact(t,x,y,z);
end

%/* ************************************************** */
function v = vely_exact(t,x,y,z)
    [u,v,w] = vel_exact(t,x,y,z);
end



