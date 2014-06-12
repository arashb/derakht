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
maxErrorPerNode = 0.001;        % Error per box
maxLevel        = 20;           % Maximum tree depth
resPerNode      = 15;           % Resolution per Node
verbose         = false;
gvfreq          = 1;
dim             = 2;
DEBUG           = false;

% MAIN SCRIPT
fconc = @func1;
fvelx = @velx;
fvely = @vely;

% CONCENTRATION TREE
c = qtree;
c.insert_function(fconc,maxErrorPerNode,maxLevel,resPerNode);
tree_init_data(c, fconc, resPerNode);

% VELOCITY TREES
u = qtree;
u.insert_function(fvelx,maxErrorPerNode,maxLevel,resPerNode);
tree_init_data(u, fvelx, resPerNode);

v = qtree;
v.insert_function(fvely,maxErrorPerNode,maxLevel,resPerNode);
tree_init_data(v, fvely, resPerNode);

% ADVECT
c_atT = advect(c, {u,u,u,u}, {v,v,v,v});

figure('Name','SEMI-LAG QUAD-TREES');
MS='MarkerSize';

subplot(3,3,2)
c.plottree;
hold on;
[cxx,cyy,cvv] = tree_griddata(c);
plot(cxx(:),cyy(:),'ro',MS,1); hold on;
axis off; axis equal;
title('c(t)');

% subplot(3,3,3)
% c_atT.plottree
% hold on;
% [ctxx,ctyy,ctvv] = tree_griddata(c_atT);
% plot(ctxx(:),ctyy(:),'ro',MS,1); hold on;
% axis off; axis equal;
% title('c(t+dt)');

subplot(3,3,5)
u.plottree;
hold on;
[uxx,uyy,uvv] = tree_griddata(u);
plot(uxx(:),uyy(:),'bo',MS,1); hold on;
axis off; axis equal;
title('u(t)');

subplot(3,3,8)
v.plottree;
hold on;
[vxx,vyy,vvv] = tree_griddata(v);
plot(vxx(:),vyy(:),'bo',MS,1); hold on;
axis off; axis equal;
title('v(t)');

figure('Name','SEMI-LAG ADVECTION')
subplot(1,2,1);
plot3(cxx,cyy,cvv,'.',MS,1);
title('c(t)');

% subplot(1,2,2);
% plot3(ctxx,ctyy,vq,'.',MS,1);
% title('c(t+dt)');
end

%/* ************************************************** */
function [c_atT] = advect(ctree, ucells, vcells)
global DEBUG;
global maxErrorPerNode;
global maxLevel;
global resPerNode;

% INIT THE C_atT
% same structure as ctree
% TODO: It will be changed later
c_atT = qtree;
c_atT.insert_function(@func1,maxErrorPerNode,maxLevel,resPerNode);

% MERGE VELOCITY TREES
um = merge(ucells);
vm = merge(vcells);

if DEBUG
    n = length(ucells);
    figure('Name','MERGED VEL. TREES');
    for i=1:n
        subplot(2,n,i)
        ucells{i}.plottree;
        axis off; axis equal;
    end
    subplot(2,n,n+1)
    um.plottree;
    axis off; axis equal;
end

% INTERPOLATE VELOCITY VALUES ON THE MERGED TREE
tree_interp(ucells{1}, um);

% CONCENTRATION VALUES
cleaves = ctree.leaves();
for lvcnt = 1:length(cleaves)
    cleaf = cleaves{lvcnt};
    %advect_sl_rk2(c, xx, yy, zz, u, v, w, t, tstep,dt, 'spline', 'spline')
end

    %/* ************************************************** */
    function mt = merge(treecells)
        mt = treecells{1};
        for counter=2:length(treecells)
            mt = qtree.merge(mt, treecells{counter});
        end
    end
end

%/* ************************************************** */
function value = func1(x,y)
xc = 0.5;
yc = 0.5;
theta = 0;
sigma = 0.05;
value = gaussian(x,y,xc,yc,theta, sigma, sigma);
end

%/* ************************************************** */
function value = velx(x,y)
t = 0;
z = 0;
xc = 0.5*ones(size(x));
yc = xc; zc = xc;
[value,v,w] = vel_rot(t,x,y,z,xc,yc,zc);
end

%/* ************************************************** */
function value = vely(x,y)
t = 0;
z = 0;
xc = 0.5*ones(size(x));
yc = xc; zc = xc;
[u,value,w] = vel_rot(t,x,y,z,xc,yc,zc);
end

%/* ************************************************** */
function value = vel(x,y)
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


