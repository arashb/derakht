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
fconc   = @conc_exact;
fvelx   = @velx_exact;
fvely   = @vely_exact;

% CONSTRUCT AND INIT THE TREES VALUES
% CONCENTRATION
c = qtree;
c.insert_function(fconc,@do_refine);
tree_data.init_data(c,fconc,resPerNode);

cdepth = c.find_depth();
cwidth = 1/2^cdepth;
dx     = cwidth/resPerNode;
cfl = 100;
om  = 1;
dt  = cfl*dx/om;
t   = [-dt 0 dt 2*dt];

% VELOCITY (TIME-DEPENDENT)
nt = length(t);
ucells = cell(1,nt);
vcells = cell(1,nt);

for tcnt = 1:nt
    utmptree = qtree;
    utmptree.insert_function(fvelx,@do_refine,t(tcnt));
    tree_data.init_data(utmptree,fvelx,resPerNode,t(tcnt));
    ucells{tcnt} = utmptree;
    
    vtmptree = qtree;
    vtmptree.insert_function(fvely,@do_refine,t(tcnt));
    tree_data.init_data(vtmptree,fvely,resPerNode,t(tcnt));
    vcells{tcnt} = vtmptree;
end

% VELOCITY (TIME-INDEPENDENT)
% u = qtree;
% u.insert_function(fvel_valx,maxErrorPerNode,maxLevel,resPerNode);
% 
% v = qtree;
% v.insert_function(fvel_valy,maxErrorPerNode,maxLevel,resPerNode);
% tree_data.init_data(u,fvel_valx,resPerNode);
% tree_data.init_data(v,fvel_valy,resPerNode);
% ucells = {u,u,u,u};
% vcells = {v,v,v,v};

% ADVECT
% INIT THE TREE FOR THE NEXT TIME STEP
% same structure as ctree
% TODO: It will be changed later
ctree_next = qtree.clone(c);

cnext = advect(c, ctree_next, ucells, vcells, t);

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

saveas(f, 'results','pdf');
    
    %/* ************************************************** */
    function val = do_refine(qtree,func,t)
        val = tree_do_refine(qtree, func, maxErrorPerNode, maxLevel, resPerNode,t);
    end
end


%/* ************************************************** */
function value = conc_exact(t,x,y,z)
xc = 0.75;
yc = 0.75;
theta = 0;
sigmax = 0.05;
sigmay = 0.05;
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
function u = velx_exact(t,x,y,z)
    [u,v,w] = vel_exact(t,x,y,z);
end

%/* ************************************************** */
function v = vely_exact(t,x,y,z)
    [u,v,w] = vel_exact(t,x,y,z);
end



