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
DEBUG           = true;
INTERP_TYPE     = 'cubic';

% MAIN SCRIPT
fconc       = @conc_tree;
fvel_valx   = @velx;
fvel_valy   = @vely;

% CONCENTRATION TREE
c = qtree;
c.insert_function(fconc,maxErrorPerNode,maxLevel,resPerNode);
tree_data.init_data(c,fconc,resPerNode);

% VELOCITY TREES
u = qtree;
u.insert_function(fvel_valx,maxErrorPerNode,maxLevel,resPerNode);
tree_data.init_data(u,fvel_valx,resPerNode);

v = qtree;
v.insert_function(fvel_valy,maxErrorPerNode,maxLevel,resPerNode);
tree_data.init_data(v,fvel_valy,resPerNode);

% ADVECT
ucells = {u,u,u,u};
vcells = {v,v,v,v};
c_atT = advect(c, ucells, vcells);

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
tree_data.plot_data(c)
colorbar;
title('c(t)');

subplot(1,2,2);
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

% INTERPOLATE VELOCITY VALUES ON THE MERGED TREE
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

% ADVECT CONCENTRATION VALUES
cdepth = ctree.find_depth();
cwidth = 1/2^cdepth;
dx = cwidth/resPerNode;

cfl = 1;
om = 1;
dt = cfl*dx/om;
t = [-2*dt -dt 0 dt];
tstep = 1;

% TODO: THIS IS GOING TO WORK ONLY FOR THE VELOCITY FIELD USED IN THIS
% EXAMPLE -> GENERALIZE THE METHOD
uc = qtree.clone(ctree);
tree_data.interp(um,uc);

vc = qtree.clone(ctree);
tree_data.interp(vm,vc);

c_leaves     = ctree.leaves();
cnext_leaves = ctree_next.leaves();
uc_leaves    = uc.leaves();
vc_leaves    = vc.leaves();

for lvcnt = 1:length(c_leaves)
    c_leaf      = c_leaves{lvcnt};
    cnext_leaf  = cnext_leaves{lvcnt};
    uc_leaf     = uc_leaves{lvcnt};
    vc_leaf     = vc_leaves{lvcnt};
        
    [xx,yy,zz,dx,dy,dz] = c_leaf.mesh(resPerNode);
    c = c_leaf.data.values;
    u = uc_leaf.data.values;
    v = vc_leaf.data.values;
    w = zeros(size(u));
    
    cnext_leaf.data.dim = c_leaf.data.dim;
    cnext_leaf.data.resolution = c_leaf.data.resolution;
    cnext_values = semilag_rk2(c,xx,yy,zz,@interp_conc, ...
        u,v,w,t,@interp_vel_precomputed,tstep,INTERP_TYPE);
    cnext_leaf.data.values(:,:,:) = cnext_values;
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
function value = conc_tree(x,y)
t=0; 
z=0;
value = conc(t,x,y,z);
end

%/* ************************************************** */
function value = conc(t,x,y,z)
xc = 0.5;
yc = 0.5;
theta = 0;
sigmax = 0.05;
sigmay = 0.09;
value = gaussian(x,y,xc,yc,theta, sigmax, sigmay);
end

%/* ************************************************** */
function [u,v,w] = vel(t,x,y,z)
z = 0;
xc = 0.5*ones(size(x));
yc = xc; zc = xc;
[u,v,w] = vel_rot(t,x,y,z,xc,yc,zc);
end

%/* ************************************************** */
function value = velx(x,y)
    [u,v,w] = vel_values(x,y);
    value = u(:,:,:,2);
end

%/* ************************************************** */
function value = vely(x,y)
    [u,v,w] = vel_values(x,y);
    value = v(:,:,:,2);
end

%/* ************************************************** */
function [u,v,w] = vel_values(x,y)
z = 0;
xc = 0.5*ones(size(x));
yc = xc; zc = xc;
dt = 0.01;
ti = 0;
tlist = linspace(ti-dt,ti+2*dt,4);
u = zeros([size(x) 1 length(tlist)]);
v = u; w = u;
for tcnt = 1:length(tlist)
    t = tlist(tcnt);
    [u(:,:,:,tcnt),v(:,:,:,tcnt),w(:,:,:,tcnt)] = vel_rot(t,x,y,z,xc,yc,zc);
end
end

%/* ************************************************** */
function value = velnorm(x,y)
t = 0;
z = 0;
value = zeros(size(x));
xc = 0.5*ones(size(x));
yc = xc; zc = xc;
[u,v,w] = vel_rot(t,x,y,z,xc,yc,zc);
for i =1:size(x,1)
    for j=1:size(y,2)
        value(i,j) = norm([u(i,j), v(i,j)]);
    end
end
end


