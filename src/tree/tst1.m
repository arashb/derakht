function tst1()
% test the creation of quadtree based on inserted function(s)
clear all; clear globals; dim=2;  % constants  and preamble
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
Z1 = f1(xx,yy); 
Z2 = f2(xx,yy);
Z3 = f3(xx,yy);
contour(xx,yy,Z1);
hold on
contour(xx,yy,Z2);
contour(xx,yy,Z3);
colorbar;
axis off;
hold on;

% create and plot the tree
o = qtree;
o.insert_function(f1,maxErrorPerNode,maxLevel,resPerNode);
o.insert_function(f2,maxErrorPerNode,maxLevel,resPerNode);
o.insert_function(f3,maxErrorPerNode,maxLevel,resPerNode);
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

    function value = func1(x,y)
        xc = 0.75;
        yc = 0.75;
        value = gaussian(x,y,xc,yc);
    end

    function value = func2(x,y)        
        xc = 0.25;
        yc = 0.25;
        value = gaussian(x,y,xc,yc);
    end

    function value = func3(x,y)
        xc = 0.5;
        yc = 0.5;
        value = gaussian(x,y,xc,yc);
    end
end

