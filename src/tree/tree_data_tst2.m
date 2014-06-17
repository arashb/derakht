function tree_data_tst2()
% test of collapssing the values of two trees with the same structure
close all; clear; clear globals; dim=2;  % constants  and preamble
addpath('./../common/');
global verbose;
global gvfreq;

% RUN PARAMETERS
maxErrorPerNode = 0.01;      % Error per box
maxLevel        = 20;          % maximum tree depth
verbose         = false;
resPerNode      = 10;
gvfreq          = 1;

% MAIN SCRIPT

% create the first tree
a = qtree;
a.insert_function(@func1,maxErrorPerNode,maxLevel,resPerNode);

% init the fist tree's data with a given function
tree_data.init_data(a,@func1,resPerNode)

% create the second tree 
b = qtree;
b.insert_function(@func1,maxErrorPerNode,maxLevel,resPerNode)

% init the fist tree's data with a given function
tree_data.init_data(b,@func2,resPerNode)

% interpolate second tree's data from the first tree
%tree_data.interp(a,b);
c = tree_data.collapse({a,b});

if verbose
   % print morton ids, all nodes
   disp(' all nodes');
   a.print_mids;
   % print morton ids, leaves only
   disp('  leaves only');
   a.print_mids(true);
end
depth=find_depth(a);
fprintf('tree depth is %d\n', depth);

% PLOTTING

subplot(3,2,1)
tree_data.plot_grid(a)
title('a grid')

subplot(3,2,2)
tree_data.plot_grid(b)
title('b grid')

subplot(3,2,3)
tree_data.plot_data(a,1)
title('a values')

subplot(3,2,4)
tree_data.plot_data(b,1)
title('b values')

subplot(3,2,5)
tree_data.plot_data(c,1)
title('c values(:,1)')

subplot(3,2,6)
tree_data.plot_data(c,2)
title('c values(:,2)')

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
        value(:,:,:,1) = u;
        value(:,:,:,2) = v;
    end
end


