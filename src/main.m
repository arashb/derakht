function main()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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
verbose         = false;
gvfreq          = 1;
dim             = 2;
DEBUG           = false;
INTERP_TYPE     = 'cubic';

% MAIN SCRIPT
fconc       = @conc_tree;
fvel_valx   = @velx_tree;
fvel_valy   = @vely_tree;

% CONSTRUCT THE TREES
% CONCENTRATION
c = qtree;
c.insert_function(fconc,maxErrorPerNode,maxLevel,resPerNode);

% VELOCITY
u = qtree;
u.insert_function(fvel_valx,maxErrorPerNode,maxLevel,resPerNode);

v = qtree;
v.insert_function(fvel_valy,maxErrorPerNode,maxLevel,resPerNode);

% INIT THE TREES VALUES
tree_data.init_data(c,fconc,resPerNode);
tree_data.init_data(u,fvel_valx,resPerNode);
tree_data.init_data(v,fvel_valy,resPerNode);

% ADVECT
ucells = {u,u,u,u};
vcells = {v,v,v,v};
c_atT = advect(c, ucells, vcells);

% PLOT THE RESULTS
figure('Name','SEMI-LAG QUAD-TREES');
subplot(3,3,2)
tree_data.plot_grid(c);
title('c(t)');

subplot(3,3,3)
tree_data.plot_grid(c_atT);
title('c(t+dt)');

subplot(3,3,5)
tree_data.plot_grid(u);
title('u(t)');

subplot(3,3,8)
tree_data.plot_grid(v);
title('v(t)');

f = figure('Name','SEMI-LAG ADVECTION');
subplot(1,2,1);
c.plottree(0.5);
tree_data.plot_data(c)
colorbar;
title('c(t)');

subplot(1,2,2);
c_atT.plottree(0.5);
tree_data.plot_data(c_atT);
colorbar;
title('c(t+dt)');

saveas(f, 'resutls','eps');
end

%/* ************************************************** */
function [ctree_next] = advect(ctree, ucells, vcells)
global DEBUG;
global resPerNode;
global INTERP_TYPE;

% INIT THE C_atT
% same structure as ctree
% TODO: It will be changed later
ctree_next = qtree.clone(ctree);

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
        tree_data.plot_data(ucells{i});
        subplot(3,n,n+i)
        tree_data.plot_data(umcells{i});
    end
    subplot(3,n,2*n+1);
    um_tree.plottree
    axis off; axis equal;
end

% PERFORM ONE STEP SEMI-LAGRANGIAN FOR EACH TREE LEAF
cdepth = ctree.find_depth();
cwidth = 1/2^cdepth;
dx     = cwidth/resPerNode;

cfl = 100;
om  = 1;
dt  = cfl*dx/om;
t   = [-dt 0 dt 2*dt];
tstep = 1;

c_leaves     = ctree.leaves();
cnext_leaves = ctree_next.leaves();
for lvcnt = 1:length(c_leaves)
    c_leaf      = c_leaves{lvcnt};
    cnext_leaf  = cnext_leaves{lvcnt};
    
    [xx,yy,zz,dx,dy,dz] = c_leaf.mesh(resPerNode);
    cnumsol = c_leaf.data.values;
    
    %fconc = @conc_exact;
    %fvel  = @vel_exact;
    fconc = @conc;    
    fvel  = @vel;
    cnext_values = semilag_rk2(xx,yy,zz,fconc,fvel,t,tstep);
    
    cnext_leaf.data.dim           = c_leaf.data.dim;
    cnext_leaf.data.resolution    = c_leaf.data.resolution;    
    cnext_leaf.data.values(:,:,:) = cnext_values;
end
    
    %/* ************************************************** */
    function ci = conc(tstep,xt,yt,zt)
        ci = interp_conc_spatial(cnumsol,xx,yy,zz,tstep,xt,yt,zt,INTERP_TYPE,@conc_out);
        
        function cq = conc_out(cq,xq,yq,zq)                 
            % OUTSIDE THE CURRENT NODE
            [xmin,xmax,ymin,ymax] = c_leaf.corners();
            out = xq<xmin | xq>xmax  | yq<ymin | yq>ymax;% | zq<0 | zq>1;
            %ce = conc_exact(0,xq,yq,zq);
            ce = tree_data.interp_points(ctree,xq,yq,zq);
            cq(out) = ce(out);
            
            % OUTSIDE THE SIMULATION DOMAIN
            out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
            cq(out) = 0;
        end
    end

    %/* ************************************************** */
    function [uq,vq,wq] = vel(tq,xq,yq,zq)
        uval = tree_data.interp_points(um,xq,yq,zq);
        vval = tree_data.interp_points(vm,xq,yq,zq);
        % TODO: THIS WORKS ONLY FOR TIME-INDEPENDENT VELOCITY FIELDS
        uq = uval(:,:,:,1);
        vq = vval(:,:,:,1);
        wq = zeros(size(uq));
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

%/* ************************************************** */
function value = conc_exact(t,x,y,z)
xc = 0.5;
yc = 0.5;
theta = 0;
sigmax = 0.05;
sigmay = 0.09;
value = gaussian(x,y,xc,yc,theta, sigmax, sigmay);
end

%/* ************************************************** */
function [u,v,w] = vel_exact(t,x,y,z)
    xc = 0.5;
    yc = 0.5;
    zc = 0.5;
    [u,v,w] = vel_rot(t,x,y,z,xc,yc,zc);
end

%/* ************************************************** */
function value = conc_tree(x,y)
t=0; 
z=0;
value = conc_exact(t,x,y,z);
end

%/* ************************************************** */
function u = velx_tree(x,y)
    [u,v,w] = vel_exact(0,x,y,0);
end

%/* ************************************************** */
function v = vely_tree(x,y)
    [u,v,w] = vel_exact(0,x,y,0);
end

%/* ************************************************** */
% function [u,v,w] = vel_values(x,y)
% z = 0;
% xc = 0.5*ones(size(x));
% yc = xc; zc = xc;
% dt = 0.01;
% ti = 0;
% tlist = linspace(ti-dt,ti+2*dt,4);
% u = zeros([size(x) 1 length(tlist)]);
% v = u; w = u;
% for tcnt = 1:length(tlist)
%     t = tlist(tcnt);
%     [u(:,:,:,tcnt),v(:,:,:,tcnt),w(:,:,:,tcnt)] = vel_rot(t,x,y,z,xc,yc,zc);
% end
% end

%/* ************************************************** */
% function value = velnorm(x,y)
% t = 0;
% z = 0;
% value = zeros(size(x));
% xc = 0.5*ones(size(x));
% yc = xc; zc = xc;
% [u,v,w] = vel_rot(t,x,y,z,xc,yc,zc);
% for i =1:size(x,1)
%     for j=1:size(y,2)
%         value(i,j) = norm([u(i,j), v(i,j)]);
%     end
% end
%end


