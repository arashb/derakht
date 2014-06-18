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

% RUN PARAMETERS
maxErrorPerNode = 0.0001;        % Error per box
maxLevel        = 20;           % Maximum tree depth
resPerNode      = 15;           % Resolution per Node
verbose         = false;
gvfreq          = 1;
dim             = 2;
DEBUG           = true;

% MAIN SCRIPT
fconc       = @func1;
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

% subplot(3,3,3)
% tree_data.plot_grid(c_atT);
% title('c(t+dt)');

subplot(3,3,5)
tree_data.plot_grid(u);
title('u(t)');

subplot(3,3,8)
tree_data.plot_grid(v);
title('v(t)');

figure('Name','SEMI-LAG ADVECTION')
subplot(1,2,1);
tree_data.plot_data(c)
title('c(t)');

% subplot(1,2,2);
%tree_data.plot_values(c_atT);
% title('c(t+dt)');
end

%/* ************************************************** */
function [c_atT] = advect(ctree, ucells, vcells)
global DEBUG;

% INIT THE C_atT
% same structure as ctree
% TODO: It will be changed later
c_atT = qtree.clone(ctree);

% MERGE VELOCITY TREES
num     = size(ucells,2);

um      = merge(ucells);
umcells = clone(um,num);

vm      = merge(vcells);
vmcells = clone(vm,num);

% INTERPOLATE VELOCITY VALUES ON THE MERGED TREE
for i =1:num
    tree_data.interp(ucells{i}, umcells{i});
    tree_data.interp(vcells{i}, vmcells{i});
end

um_col = tree_data.collapse(umcells);
vm_col = tree_data.collapse(vmcells);

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
    um.plottree
    axis off; axis equal;
end

% CONCENTRATION VALUES
cleaves     = ctree.leaves();
cTleaves    = c_atT.leaves();
for lvcnt = 1:length(cleaves)
    cleaf  = cleaves{lvcnt};
    cTleaf = cTleaves{lvcnt};
    % = advect_sl_rk2(c, xx, yy, zz, u, v, w, t, tstep,dt, 'spline', 'spline')
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
function value = func1(x,y)
xc = 0.5;
yc = 0.5;
theta = 0;
sigmax = 0.05;
sigmay = 0.09;
value = gaussian(x,y,xc,yc,theta, sigmax, sigmay);
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


