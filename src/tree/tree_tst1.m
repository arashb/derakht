function tree_tst1()
% test the creation of quadtree based on inserted function(s)
clear all; clear globals; dim=2;  % constants  and preamble
addpath('../common/');
global verbose;

% RUN PARAMETERS
maxErrorPerNode = 0.0001;      % Error per box
maxLevel        = 20;          % maximum tree depth
resPerNode      = 10;          % Resolution per box
verbose         = false;

% MAIN SCRIPT
f1 = @func1;
f2 = @func2;
f3 = @func3;

% plot the function
res = 1000;
x = linspace(0,1,res);
y = linspace(0,1,res);
[xx,yy] = meshgrid(x,y);
Z1 = f1(0,xx,yy,0); 
Z2 = f2(0,xx,yy,0);
Z3 = f3(0,xx,yy,0);
contour(xx,yy,Z1);
hold on
contour(xx,yy,Z2);
contour(xx,yy,Z3);
colorbar;
axis off;
hold on;

% create and plot the tree
o = qtree;
o.insert_function(f1,@do_refine);
o.insert_function(f2,@do_refine);
o.insert_function(f3,@do_refine);
o.plottree;
axis off; hold on;

if verbose
   % print morton ids, all nodes
   disp(' all nodes');
   o.print_mids;
   % print morton ids, leaves only
   disp('  leaves only');
   o.print_mids(true);
end
depth=find_depth(o);
fprintf('tree depth is %d\n', depth);

    function value = func1(t,x,y,z)
        xc = 0.75;
        yc = 0.75;
        value = gaussian(x,y,xc,yc);
    end

    function value = func2(t,x,y,z)        
        xc = 0.25;
        yc = 0.25;
        value = gaussian(x,y,xc,yc);
    end

    function value = func3(t,x,y,z)
        xc = 0.5;
        yc = 0.5;
        value = gaussian(x,y,xc,yc);
    end

    function val = do_refine(qtree,func,t)
        val = tree_do_refine(qtree, func, maxErrorPerNode, maxLevel, resPerNode,t);
    end
end

