function tst4()
% test the creation of quadtree based on inserted function(s)
clear all; clear globals; dim=2;  % constants  and preamble
global verbose;
global gvfreq;

% RUN PARAMETERS
maxErrorPerNode = 0.00001;      % Error per box
maxLevel        = 20;          % maximum tree depth
verbose         = false;
gvfreq           = 1;

% MAIN SCRIPT

% plot the function
res = 1000;
x = linspace(0,1,res);
y = linspace(0,1,res);
[xx,yy] = meshgrid(x,y);
Z1 = func3(xx,yy); 

contour(xx,yy,Z1);
colorbar;
axis off;
hold on;

% create and plot the tree
o = qtree;
o.insert_function(@func3,maxErrorPerNode,maxLevel);
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
end

